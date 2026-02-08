#!/usr/bin/env python3
"""
FastQC 결과 자동 평가 스크립트

FastQC 분석 결과를 파싱하고 FASTQC_GUIDE.md에 명시된 기준에 따라
각 샘플의 품질을 자동으로 판단합니다.

RNA-seq 특성을 고려하여 정상적인 FAIL 패턴을 허용하고,
실제 문제가 있는 샘플만 플래그합니다.
"""

import sys
import os
import re
import json
import argparse
from pathlib import Path
from typing import Dict, List, Tuple
from collections import defaultdict


class FastQCEvaluator:
    """FastQC 결과를 평가하는 클래스"""
    
    # RNA-seq에서 정상적으로 FAIL이 나올 수 있는 모듈들
    RNA_SEQ_ACCEPTABLE_FAILS = {
        'Per base sequence content',
        'Sequence Duplication Levels',
        'Overrepresented sequences',
        'Adapter Content'  # trimming 전에만
    }
    
    def __init__(self, config: Dict):
        """
        Args:
            config: QC 평가 기준을 담은 설정 딕셔너리
        """
        self.config = config
        self.results = []
        
    def parse_fastqc_data(self, fastqc_data_path: str) -> Dict:
        """
        FastQC의 fastqc_data.txt 파일을 파싱
        
        Args:
            fastqc_data_path: fastqc_data.txt 파일 경로
            
        Returns:
            파싱된 데이터 딕셔너리
        """
        data = {
            'basic_stats': {},
            'per_base_quality': [],
            'per_sequence_quality': [],
            'per_base_n_content': [],
            'sequence_length': [],
            'adapter_content': [],
            'modules': {}
        }
        
        current_module = None
        in_data_section = False
        
        with open(fastqc_data_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                
                if line.startswith('>>') and not line.startswith('>>END_MODULE'):
                    # 새로운 모듈 시작
                    parts = line[2:].split('\t')
                    module_name = parts[0]
                    module_status = parts[1] if len(parts) > 1 else 'unknown'
                    current_module = module_name
                    data['modules'][module_name] = module_status
                    in_data_section = True
                    continue
                    
                if line.startswith('>>END_MODULE'):
                    current_module = None
                    in_data_section = False
                    continue
                
                if line.startswith('#') or not line:
                    continue
                
                if not in_data_section:
                    continue
                
                # 각 모듈별 데이터 파싱
                if current_module == 'Basic Statistics':
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        data['basic_stats'][parts[0]] = parts[1]
                        
                elif current_module == 'Per base sequence quality':
                    parts = line.split('\t')
                    if len(parts) >= 7:
                        data['per_base_quality'].append({
                            'base': parts[0],
                            'mean': float(parts[1]),
                            'median': float(parts[2]),
                            'lower_quartile': float(parts[3]),
                            'upper_quartile': float(parts[4]),
                            'percentile_10': float(parts[5]),
                            'percentile_90': float(parts[6])
                        })
                        
                elif current_module == 'Per sequence quality scores':
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        data['per_sequence_quality'].append({
                            'quality': int(parts[0]),
                            'count': float(parts[1])
                        })
                        
                elif current_module == 'Per base N content':
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        data['per_base_n_content'].append({
                            'base': parts[0],
                            'n_count': float(parts[1])
                        })
                        
                elif current_module == 'Adapter Content':
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        # 모든 어댑터 컬럼의 최대값 추출
                        max_adapter = max([float(x) for x in parts[1:] if x])
                        data['adapter_content'].append({
                            'position': parts[0],
                            'max_percentage': max_adapter
                        })
        
        return data
    
    def evaluate_sample(self, sample_name: str, fastqc_zip_path: str, 
                       is_trimmed: bool = False) -> Dict:
        """
        개별 샘플의 FastQC 결과를 평가
        
        Args:
            sample_name: 샘플 이름
            fastqc_zip_path: FastQC zip 파일 경로
            is_trimmed: trimming 후 데이터인지 여부
            
        Returns:
            평가 결과 딕셔너리
        """
        # zip 파일 압축 해제 경로 찾기
        zip_path = Path(fastqc_zip_path)
        data_dir = zip_path.parent / zip_path.stem
        fastqc_data_path = data_dir / 'fastqc_data.txt'
        
        if not fastqc_data_path.exists():
            return {
                'sample': sample_name,
                'status': 'ERROR',
                'message': f'fastqc_data.txt not found: {fastqc_data_path}',
                'issues': [],
                'warnings': []
            }
        
        # 데이터 파싱
        data = self.parse_fastqc_data(str(fastqc_data_path))
        
        # 평가 수행
        issues = []
        warnings = []
        
        # 1. Basic Statistics 확인
        basic_issues = self._check_basic_stats(data['basic_stats'])
        issues.extend(basic_issues)
        
        # 2. Per base sequence quality 확인
        quality_result = self._check_per_base_quality(data['per_base_quality'])
        if quality_result['severity'] == 'critical':
            issues.append(quality_result['message'])
        elif quality_result['severity'] == 'warning':
            warnings.append(quality_result['message'])
        
        # 3. Per base N content 확인
        n_content_result = self._check_n_content(data['per_base_n_content'])
        if n_content_result['severity'] == 'critical':
            issues.append(n_content_result['message'])
        elif n_content_result['severity'] == 'warning':
            warnings.append(n_content_result['message'])
        
        # 4. Adapter Content 확인 (trimming 후라면 필수)
        if is_trimmed:
            adapter_result = self._check_adapter_content(data['adapter_content'], is_trimmed)
            if adapter_result['severity'] == 'critical':
                issues.append(adapter_result['message'])
            elif adapter_result['severity'] == 'warning':
                warnings.append(adapter_result['message'])
        
        # 5. 모듈 상태 확인 (RNA-seq 정상 패턴 제외)
        module_issues = self._check_modules(data['modules'], is_trimmed)
        warnings.extend(module_issues)
        
        # 6. Per sequence quality 확인
        seq_quality_result = self._check_per_sequence_quality(data['per_sequence_quality'])
        if seq_quality_result['severity'] == 'warning':
            warnings.append(seq_quality_result['message'])
        
        # 최종 판정
        if issues:
            status = 'FAIL'
            recommendation = 'REVIEW REQUIRED - 심각한 문제 발견. 재시퀀싱 또는 제외 고려'
        elif warnings:
            status = 'WARN'
            recommendation = 'CHECK RECOMMENDED - 경미한 문제 있음. 리포트 확인 권장'
        else:
            status = 'PASS'
            recommendation = 'GOOD - 분석 사용 가능'
        
        return {
            'sample': sample_name,
            'status': status,
            'recommendation': recommendation,
            'basic_stats': data['basic_stats'],
            'issues': issues,
            'warnings': warnings,
            'modules': data['modules']
        }
    
    def _check_basic_stats(self, stats: Dict) -> List[str]:
        """기본 통계 확인"""
        issues = []
        
        # Total Sequences 확인
        if 'Total Sequences' in stats:
            total_seq = int(stats['Total Sequences'].replace(',', ''))
            min_sequences = self.config.get('min_total_sequences', 1000000)
            
            if total_seq < min_sequences:
                issues.append(
                    f'Total Sequences too low: {total_seq:,} (minimum: {min_sequences:,})'
                )
        
        # GC content 확인
        if '%GC' in stats:
            gc_content = int(stats['%GC'])
            min_gc = self.config.get('min_gc_content', 30)
            max_gc = self.config.get('max_gc_content', 70)
            
            if gc_content < min_gc or gc_content > max_gc:
                issues.append(
                    f'GC content out of range: {gc_content}% (expected: {min_gc}-{max_gc}%)'
                )
        
        return issues
    
    def _check_per_base_quality(self, quality_data: List[Dict]) -> Dict:
        """염기별 품질 점수 확인"""
        if not quality_data:
            return {'severity': 'ok', 'message': ''}
        
        min_median_quality = self.config.get('min_median_quality', 28)
        min_lower_quartile = self.config.get('min_lower_quartile', 20)
        
        critical_positions = []
        warning_positions = []
        
        for item in quality_data:
            position = item['base']
            median = item['median']
            lower_q = item['lower_quartile']
            
            if median < 20:
                critical_positions.append(f"{position} (Q{median:.1f})")
            elif median < min_median_quality:
                warning_positions.append(f"{position} (Q{median:.1f})")
            elif lower_q < 10:
                critical_positions.append(f"{position} (LQ{lower_q:.1f})")
        
        if critical_positions:
            return {
                'severity': 'critical',
                'message': f'Critical: Low quality at positions: {", ".join(critical_positions[:5])}'
            }
        elif warning_positions:
            return {
                'severity': 'warning',
                'message': f'Warning: Quality below Q{min_median_quality} at {len(warning_positions)} positions'
            }
        
        return {'severity': 'ok', 'message': ''}
    
    def _check_per_sequence_quality(self, quality_data: List[Dict]) -> Dict:
        """시퀀스별 품질 점수 분포 확인"""
        if not quality_data:
            return {'severity': 'ok', 'message': ''}
        
        # 전체 리드 수 계산
        total_reads = sum(item['count'] for item in quality_data)
        
        # Q30 이상 리드 비율 계산
        high_quality_reads = sum(
            item['count'] for item in quality_data 
            if item['quality'] >= 30
        )
        
        q30_percentage = (high_quality_reads / total_reads * 100) if total_reads > 0 else 0
        min_q30_percentage = self.config.get('min_q30_percentage', 75)
        
        if q30_percentage < min_q30_percentage:
            return {
                'severity': 'warning',
                'message': f'Low Q30+ reads: {q30_percentage:.1f}% (expected: >{min_q30_percentage}%)'
            }
        
        return {'severity': 'ok', 'message': ''}
    
    def _check_n_content(self, n_content_data: List[Dict]) -> Dict:
        """N 염기 함량 확인"""
        if not n_content_data:
            return {'severity': 'ok', 'message': ''}
        
        max_n_content = self.config.get('max_n_content', 5.0)
        critical_n_content = self.config.get('critical_n_content', 10.0)
        
        max_n = max(item['n_count'] for item in n_content_data)
        
        if max_n > critical_n_content:
            return {
                'severity': 'critical',
                'message': f'Critical: High N content: {max_n:.1f}% (max allowed: {critical_n_content}%)'
            }
        elif max_n > max_n_content:
            return {
                'severity': 'warning',
                'message': f'Warning: N content: {max_n:.1f}% (recommended: <{max_n_content}%)'
            }
        
        return {'severity': 'ok', 'message': ''}
    
    def _check_adapter_content(self, adapter_data: List[Dict], is_trimmed: bool) -> Dict:
        """어댑터 함량 확인"""
        if not adapter_data:
            return {'severity': 'ok', 'message': ''}
        
        max_adapter = max(item['max_percentage'] for item in adapter_data)
        
        if is_trimmed:
            # Trimming 후에는 엄격한 기준
            max_allowed = self.config.get('max_adapter_trimmed', 1.0)
            if max_adapter > max_allowed:
                return {
                    'severity': 'critical',
                    'message': f'Adapter content after trimming: {max_adapter:.2f}% (should be <{max_allowed}%)'
                }
        else:
            # Trimming 전에는 경고만
            warn_threshold = self.config.get('warn_adapter_raw', 10.0)
            if max_adapter > warn_threshold:
                return {
                    'severity': 'warning',
                    'message': f'Adapter content: {max_adapter:.2f}% - Trimming recommended'
                }
        
        return {'severity': 'ok', 'message': ''}
    
    def _check_modules(self, modules: Dict, is_trimmed: bool) -> List[str]:
        """
        각 모듈의 PASS/WARN/FAIL 상태 확인
        RNA-seq에서 정상적인 FAIL은 제외
        """
        warnings = []
        
        for module_name, status in modules.items():
            # RNA-seq에서 정상적으로 FAIL이 나는 모듈은 건너뜀
            if module_name in self.RNA_SEQ_ACCEPTABLE_FAILS:
                # Adapter Content는 trimming 후에는 체크해야 함
                if module_name == 'Adapter Content' and is_trimmed:
                    if status == 'fail':
                        warnings.append(
                            f'Module FAIL: {module_name} (should pass after trimming)'
                        )
                continue
            
            # 나머지 모듈에서 FAIL은 경고
            if status == 'fail':
                # Per tile sequence quality는 흔하므로 경고만
                if module_name == 'Per tile sequence quality':
                    continue  # 무시
                else:
                    warnings.append(f'Module FAIL: {module_name}')
        
        return warnings


def generate_summary_report(results: List[Dict], output_path: str):
    """
    모든 샘플의 평가 결과를 요약한 리포트 생성
    
    Args:
        results: 각 샘플의 평가 결과 리스트
        output_path: 출력 파일 경로
    """
    # 상태별 샘플 분류
    status_groups = defaultdict(list)
    for result in results:
        status_groups[result['status']].append(result)
    
    # 리포트 생성
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write("# FastQC 자동 평가 리포트\n\n")
        f.write("=" * 80 + "\n\n")
        
        # 요약
        f.write("## 📊 요약\n\n")
        f.write(f"- 총 샘플 수: {len(results)}\n")
        f.write(f"- ✅ PASS: {len(status_groups['PASS'])} 샘플\n")
        f.write(f"- ⚠️ WARN: {len(status_groups['WARN'])} 샘플\n")
        f.write(f"- ❌ FAIL: {len(status_groups['FAIL'])} 샘플\n")
        f.write(f"- 🔴 ERROR: {len(status_groups['ERROR'])} 샘플\n\n")
        
        f.write("=" * 80 + "\n\n")
        
        # FAIL 샘플 (심각한 문제)
        if status_groups['FAIL']:
            f.write("## ❌ FAIL - 즉시 확인 필요\n\n")
            f.write("**심각한 품질 문제가 발견되었습니다. 재시퀀싱 또는 분석 제외를 고려하세요.**\n\n")
            
            for result in sorted(status_groups['FAIL'], key=lambda x: x['sample']):
                f.write(f"### {result['sample']}\n\n")
                f.write(f"**권장사항**: {result['recommendation']}\n\n")
                
                if result['issues']:
                    f.write("**문제점**:\n")
                    for issue in result['issues']:
                        f.write(f"- 🔴 {issue}\n")
                    f.write("\n")
                
                if result['warnings']:
                    f.write("**경고**:\n")
                    for warning in result['warnings']:
                        f.write(f"- ⚠️ {warning}\n")
                    f.write("\n")
                
                if 'basic_stats' in result:
                    f.write("**기본 통계**:\n")
                    f.write(f"- Total Sequences: {result['basic_stats'].get('Total Sequences', 'N/A')}\n")
                    f.write(f"- GC Content: {result['basic_stats'].get('%GC', 'N/A')}%\n")
                    f.write(f"- Sequence Length: {result['basic_stats'].get('Sequence length', 'N/A')}\n\n")
                
                f.write("-" * 80 + "\n\n")
        
        # WARN 샘플 (경미한 문제)
        if status_groups['WARN']:
            f.write("## ⚠️ WARN - 리포트 확인 권장\n\n")
            f.write("**경미한 문제가 있습니다. FastQC HTML 리포트를 확인하세요.**\n\n")
            
            for result in sorted(status_groups['WARN'], key=lambda x: x['sample']):
                f.write(f"### {result['sample']}\n\n")
                
                if result['warnings']:
                    f.write("**경고**:\n")
                    for warning in result['warnings']:
                        f.write(f"- ⚠️ {warning}\n")
                    f.write("\n")
                
                f.write("-" * 80 + "\n\n")
        
        # PASS 샘플 (정상)
        if status_groups['PASS']:
            f.write("## ✅ PASS - 분석 사용 가능\n\n")
            f.write("**품질 기준을 통과했습니다. 다음 분석 단계로 진행 가능합니다.**\n\n")
            
            pass_samples = sorted([r['sample'] for r in status_groups['PASS']])
            for i in range(0, len(pass_samples), 5):
                f.write("- " + ", ".join(pass_samples[i:i+5]) + "\n")
            f.write("\n")
        
        # ERROR 샘플
        if status_groups['ERROR']:
            f.write("## 🔴 ERROR - 평가 실패\n\n")
            for result in sorted(status_groups['ERROR'], key=lambda x: x['sample']):
                f.write(f"### {result['sample']}\n\n")
                f.write(f"**에러 메시지**: {result.get('message', 'Unknown error')}\n\n")
                f.write("-" * 80 + "\n\n")
        
        # 추가 안내
        f.write("=" * 80 + "\n\n")
        f.write("## 📝 참고사항\n\n")
        f.write("### RNA-seq에서 정상적인 FAIL 패턴\n\n")
        f.write("다음 모듈들은 RNA-seq 데이터의 특성상 FAIL이 나와도 정상입니다:\n\n")
        f.write("- **Per base sequence content**: 랜덤 헥사머 프라이밍 편향 (처음 10-15bp)\n")
        f.write("- **Sequence Duplication Levels**: 고발현 유전자 (Actb, Gapdh, rRNA)\n")
        f.write("- **Overrepresented sequences**: rRNA, 미토콘드리아, 하우스키핑 유전자\n")
        f.write("- **Adapter Content** (trimming 전): Insert < Read length\n\n")
        
        f.write("### 다음 단계\n\n")
        f.write("1. **FAIL 샘플**: FastQC HTML 리포트를 열어 상세 확인\n")
        f.write("2. **WARN 샘플**: 경고 내용 검토 후 진행 여부 결정\n")
        f.write("3. **PASS 샘플**: 다음 분석 단계 진행\n\n")
        
        f.write("자세한 해석 가이드는 `FASTQC_GUIDE.md`를 참조하세요.\n\n")


def main():
    parser = argparse.ArgumentParser(
        description='FastQC 결과를 자동으로 평가하여 품질 문제를 식별합니다.'
    )
    parser.add_argument(
        'qc_dir',
        help='FastQC 결과가 있는 디렉토리 (예: results/qc/)'
    )
    parser.add_argument(
        '-o', '--output',
        default='results/qc/fastqc_evaluation.txt',
        help='평가 리포트 출력 경로 (기본값: results/qc/fastqc_evaluation.txt)'
    )
    parser.add_argument(
        '--json',
        help='JSON 형식으로 결과 저장할 경로 (선택사항)'
    )
    parser.add_argument(
        '--config',
        help='QC 기준 설정 파일 (JSON, 선택사항)'
    )
    parser.add_argument(
        '--trimmed',
        action='store_true',
        help='Trimming 후 데이터 평가 (어댑터 함량 엄격히 체크)'
    )
    
    args = parser.parse_args()
    
    # 설정 로드
    config = {
        'min_total_sequences': 1000000,
        'min_gc_content': 30,
        'max_gc_content': 70,
        'min_median_quality': 28,
        'min_lower_quartile': 20,
        'min_q30_percentage': 75,
        'max_n_content': 5.0,
        'critical_n_content': 10.0,
        'max_adapter_trimmed': 1.0,
        'warn_adapter_raw': 10.0
    }
    
    if args.config and os.path.exists(args.config):
        with open(args.config, 'r') as f:
            custom_config = json.load(f)
            config.update(custom_config)
    
    # Evaluator 생성
    evaluator = FastQCEvaluator(config)
    
    # QC 디렉토리에서 FastQC 결과 찾기
    qc_dir = Path(args.qc_dir)
    if not qc_dir.exists():
        print(f"ERROR: QC directory not found: {qc_dir}", file=sys.stderr)
        sys.exit(1)
    
    # FastQC zip 파일 찾기
    fastqc_files = list(qc_dir.glob('*_fastqc.zip'))
    
    if not fastqc_files:
        print(f"ERROR: No FastQC results found in {qc_dir}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Found {len(fastqc_files)} FastQC results")
    print(f"Evaluating with RNA-seq specific criteria...")
    print(f"Trimmed data: {args.trimmed}")
    print()
    
    # 각 샘플 평가
    results = []
    for fastqc_zip in sorted(fastqc_files):
        sample_name = fastqc_zip.stem.replace('_fastqc', '')
        print(f"Evaluating: {sample_name}...", end=' ')
        
        result = evaluator.evaluate_sample(
            sample_name,
            str(fastqc_zip),
            is_trimmed=args.trimmed
        )
        results.append(result)
        
        # 상태에 따라 색상 출력 (터미널 지원시)
        status_symbol = {
            'PASS': '✅',
            'WARN': '⚠️',
            'FAIL': '❌',
            'ERROR': '🔴'
        }.get(result['status'], '?')
        
        print(f"{status_symbol} {result['status']}")
    
    print()
    print("=" * 80)
    print(f"Evaluation complete!")
    print(f"  PASS: {sum(1 for r in results if r['status'] == 'PASS')}")
    print(f"  WARN: {sum(1 for r in results if r['status'] == 'WARN')}")
    print(f"  FAIL: {sum(1 for r in results if r['status'] == 'FAIL')}")
    print(f"  ERROR: {sum(1 for r in results if r['status'] == 'ERROR')}")
    print("=" * 80)
    print()
    
    # 리포트 생성
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    generate_summary_report(results, args.output)
    print(f"Summary report saved to: {args.output}")
    
    # JSON 저장 (선택사항)
    if args.json:
        with open(args.json, 'w', encoding='utf-8') as f:
            json.dump(results, f, indent=2, ensure_ascii=False)
        print(f"JSON results saved to: {args.json}")
    
    # FAIL이 있으면 exit code 1 반환 (선택사항으로 사용 가능)
    # if any(r['status'] == 'FAIL' for r in results):
    #     sys.exit(1)
    
    sys.exit(0)


if __name__ == '__main__':
    main()
