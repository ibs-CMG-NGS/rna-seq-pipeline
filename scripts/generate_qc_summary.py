#!/usr/bin/env python3
"""
QC Summary Generator

STAR alignment 로그와 featureCounts summary를 파싱하여
표준화된 qc_summary.json 파일을 생성합니다.
"""

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional


class QCSummaryGenerator:
    """QC 요약 생성기"""
    
    def __init__(self, sample_id: str, pipeline_type: str = "rna-seq"):
        """
        Args:
            sample_id: 샘플 ID
            pipeline_type: 파이프라인 타입
        """
        self.sample_id = sample_id
        self.pipeline_type = pipeline_type
        
    def parse_star_log(self, log_path: str) -> Dict:
        """STAR Log.final.out 파일 파싱"""
        metrics = {}
        
        with open(log_path, 'r') as f:
            for line in f:
                line = line.strip()
                if '|' in line:
                    parts = [p.strip() for p in line.split('|')]
                    if len(parts) == 2:
                        key, value = parts
                        
                        # 숫자 값 추출
                        if 'Number of input reads' in key:
                            metrics['total_reads'] = int(value)
                        elif 'Uniquely mapped reads number' in key:
                            metrics['uniquely_mapped_reads'] = int(value)
                        elif 'Uniquely mapped reads %' in key:
                            metrics['uniquely_mapped_pct'] = float(value.rstrip('%'))
                        elif 'Number of reads mapped to multiple loci' in key:
                            metrics['multi_mapped_reads'] = int(value)
                        elif '% of reads mapped to multiple loci' in key:
                            metrics['multi_mapped_pct'] = float(value.rstrip('%'))
                        elif 'Number of reads unmapped: too short' in key:
                            metrics['unmapped_too_short'] = int(value)
                        elif '% of reads unmapped: too short' in key:
                            metrics['unmapped_too_short_pct'] = float(value.rstrip('%'))
                        elif 'Average input read length' in key:
                            metrics['read_length'] = int(value)
                        elif 'Mismatch rate per base, %' in key:
                            metrics['mismatch_rate'] = float(value.rstrip('%'))
        
        # Total mapped 계산
        if 'uniquely_mapped_reads' in metrics and 'multi_mapped_reads' in metrics:
            metrics['total_mapped_reads'] = metrics['uniquely_mapped_reads'] + metrics['multi_mapped_reads']
            if 'total_reads' in metrics and metrics['total_reads'] > 0:
                metrics['total_mapped_pct'] = (metrics['total_mapped_reads'] / metrics['total_reads']) * 100
        
        return metrics
    
    def parse_featurecounts_summary(self, summary_path: str, sample_id: str) -> Dict:
        """featureCounts summary 파일 파싱"""
        metrics = {}
        
        with open(summary_path, 'r') as f:
            lines = f.readlines()
        
        if len(lines) < 2:
            return metrics
        
        # 헤더에서 샘플 컬럼 찾기
        header = lines[0].strip().split('\t')
        sample_col = None
        
        for i, col in enumerate(header):
            if sample_id in col:
                sample_col = i
                break
        
        if sample_col is None:
            print(f"Warning: Sample {sample_id} not found in featureCounts summary", file=sys.stderr)
            return metrics
        
        # 데이터 파싱
        for line in lines[1:]:
            parts = line.strip().split('\t')
            if len(parts) <= sample_col:
                continue
            
            status = parts[0]
            count = int(parts[sample_col])
            
            if status == 'Assigned':
                metrics['assigned_reads'] = count
            elif status == 'Unassigned_MultiMapping':
                metrics['unassigned_multimapping'] = count
            elif status == 'Unassigned_NoFeatures':
                metrics['unassigned_nofeatures'] = count
            elif status == 'Unassigned_Ambiguity':
                metrics['unassigned_ambiguity'] = count
        
        # Assignment rate 계산
        if 'assigned_reads' in metrics:
            total = sum([
                metrics.get('assigned_reads', 0),
                metrics.get('unassigned_multimapping', 0),
                metrics.get('unassigned_nofeatures', 0),
                metrics.get('unassigned_ambiguity', 0)
            ])
            
            if total > 0:
                metrics['assignment_rate'] = (metrics['assigned_reads'] / total) * 100
            else:
                metrics['assignment_rate'] = 0.0
        
        return metrics
    
    def determine_status(self, star_metrics: Dict, fc_metrics: Dict) -> str:
        """QC 상태 판정"""
        issues = []
        
        # Mapping rate 확인
        uniquely_mapped = star_metrics.get('uniquely_mapped_pct', 0)
        if uniquely_mapped < 70:
            issues.append(f"Low uniquely mapped rate: {uniquely_mapped:.1f}%")
        
        # Assignment rate 확인
        assignment_rate = fc_metrics.get('assignment_rate', 0)
        if assignment_rate < 50:
            issues.append(f"Low assignment rate: {assignment_rate:.1f}%")
        
        # 판정
        if uniquely_mapped < 50 or assignment_rate < 30:
            return "FAIL"
        elif issues:
            return "WARN"
        else:
            return "PASS"
    
    def get_recommendations(self, status: str, issues: List[str]) -> List[str]:
        """권장사항 생성"""
        if status == "PASS":
            return ["Sample quality is sufficient for downstream analysis"]
        elif status == "WARN":
            return [
                "Review QC report for detailed metrics",
                "Consider sample in analysis with caution",
                *issues
            ]
        else:  # FAIL
            return [
                "Sample quality is below acceptable thresholds",
                "Consider excluding from analysis",
                "May require resequencing",
                *issues
            ]
    
    def generate_summary(self, 
                        star_log: str, 
                        featurecounts_summary: str,
                        fastqc_data: Optional[str] = None,
                        output_path: Optional[str] = None) -> Dict:
        """QC summary 생성"""
        
        # STAR 메트릭 파싱
        star_metrics = self.parse_star_log(star_log)
        
        # featureCounts 메트릭 파싱
        fc_metrics = self.parse_featurecounts_summary(featurecounts_summary, self.sample_id)
        
        # 상태 판정
        status = self.determine_status(star_metrics, fc_metrics)
        
        # Issues 수집
        issues = []
        if star_metrics.get('uniquely_mapped_pct', 0) < 70:
            issues.append({
                "severity": "warning" if star_metrics.get('uniquely_mapped_pct', 0) >= 50 else "critical",
                "category": "alignment",
                "message": f"Low uniquely mapped rate: {star_metrics.get('uniquely_mapped_pct', 0):.1f}%"
            })
        
        if fc_metrics.get('assignment_rate', 0) < 50:
            issues.append({
                "severity": "warning" if fc_metrics.get('assignment_rate', 0) >= 30 else "critical",
                "category": "quantification",
                "message": f"Low assignment rate: {fc_metrics.get('assignment_rate', 0):.1f}%"
            })
        
        # Summary 구성
        summary = {
            "sample_id": self.sample_id,
            "pipeline_type": self.pipeline_type,
            "qc_version": "1.0",
            "overall_status": status,
            
            "sequencing": {
                "total_reads": star_metrics.get('total_reads', 0),
                "read_length": star_metrics.get('read_length', 0),
            },
            
            "alignment": {
                "total_mapped_reads": star_metrics.get('total_mapped_reads', 0),
                "total_mapped_pct": star_metrics.get('total_mapped_pct', 0.0),
                "uniquely_mapped_reads": star_metrics.get('uniquely_mapped_reads', 0),
                "uniquely_mapped_pct": star_metrics.get('uniquely_mapped_pct', 0.0),
                "multi_mapped_reads": star_metrics.get('multi_mapped_reads', 0),
                "multi_mapped_pct": star_metrics.get('multi_mapped_pct', 0.0),
                "unmapped_too_short": star_metrics.get('unmapped_too_short', 0),
                "unmapped_too_short_pct": star_metrics.get('unmapped_too_short_pct', 0.0),
                "mismatch_rate": star_metrics.get('mismatch_rate', 0.0)
            },
            
            "quantification": {
                "assigned_reads": fc_metrics.get('assigned_reads', 0),
                "assignment_rate": fc_metrics.get('assignment_rate', 0.0),
                "unassigned_multimapping": fc_metrics.get('unassigned_multimapping', 0),
                "unassigned_nofeatures": fc_metrics.get('unassigned_nofeatures', 0),
                "unassigned_ambiguity": fc_metrics.get('unassigned_ambiguity', 0)
            },
            
            "issues": issues,
            
            "recommendations": self.get_recommendations(status, [i['message'] for i in issues])
        }
        
        # 파일로 저장
        if output_path:
            output_file = Path(output_path)
            output_file.parent.mkdir(parents=True, exist_ok=True)
            
            with open(output_file, 'w', encoding='utf-8') as f:
                json.dump(summary, f, indent=2, ensure_ascii=False)
        
        return summary


def main():
    parser = argparse.ArgumentParser(
        description='Generate qc_summary.json from pipeline logs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python generate_qc_summary.py \\
      --sample-id h_RNA_Cont_1 \\
      --star-log aligned/h_RNA_Cont_1/Log.final.out \\
      --featurecounts counts/counts_matrix.txt.summary \\
      -o qc/qc_summary.json
        """
    )
    
    parser.add_argument(
        '--sample-id',
        required=True,
        help='Sample ID'
    )
    parser.add_argument(
        '--star-log',
        required=True,
        help='STAR Log.final.out file path'
    )
    parser.add_argument(
        '--featurecounts',
        required=True,
        help='featureCounts summary file path'
    )
    parser.add_argument(
        '--pipeline-type',
        default='rna-seq',
        help='Pipeline type (default: rna-seq)'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output qc_summary.json path'
    )
    
    args = parser.parse_args()
    
    try:
        generator = QCSummaryGenerator(args.sample_id, args.pipeline_type)
        
        print(f"\n{'='*70}")
        print(f"Generating QC Summary")
        print(f"  Sample: {args.sample_id}")
        print(f"  STAR log: {args.star_log}")
        print(f"  featureCounts: {args.featurecounts}")
        print(f"{'='*70}\n")
        
        summary = generator.generate_summary(
            args.star_log,
            args.featurecounts,
            output_path=args.output
        )
        
        print(f"✅ QC summary generated: {args.output}")
        print(f"\nStatus: {summary['overall_status']}")
        print(f"  Uniquely mapped: {summary['alignment']['uniquely_mapped_pct']:.2f}%")
        print(f"  Assignment rate: {summary['quantification']['assignment_rate']:.2f}%")
        
        if summary['issues']:
            print(f"\nIssues found: {len(summary['issues'])}")
            for issue in summary['issues']:
                print(f"  {issue['severity'].upper()}: {issue['message']}")
        else:
            print(f"\n✅ No issues found")
        
        print()
        
        return 0
        
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
