#!/usr/bin/env python3
"""
결과 폴더 구조 검증 스크립트
사용법: python scripts/verify_results.py
"""

import os
import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple

class Colors:
    """터미널 색상 코드"""
    GREEN = '\033[0;32m'
    RED = '\033[0;31m'
    YELLOW = '\033[1;33m'
    BLUE = '\033[0;34m'
    NC = '\033[0m'  # No Color
    BOLD = '\033[1m'

def print_header(text: str):
    """헤더 출력"""
    print(f"\n{Colors.BOLD}{'='*60}{Colors.NC}")
    print(f"{Colors.BOLD}{text}{Colors.NC}")
    print(f"{Colors.BOLD}{'='*60}{Colors.NC}\n")

def print_section(text: str):
    """섹션 헤더 출력"""
    print(f"\n{Colors.BLUE}{text}{Colors.NC}")
    print(f"{'-'*60}")

def check_file(path: Path, name: str, critical: bool = True) -> bool:
    """파일 존재 및 크기 확인"""
    if path.exists():
        size = path.stat().st_size
        size_str = format_size(size)
        print(f"{Colors.GREEN}✓{Colors.NC} {name}: {path} ({size_str})")
        return True
    else:
        if critical:
            print(f"{Colors.RED}✗{Colors.NC} {name}: {path} (NOT FOUND)")
        else:
            print(f"{Colors.YELLOW}⚠{Colors.NC} {name}: {path} (OPTIONAL - NOT FOUND)")
        return False

def check_dir(path: Path, name: str, critical: bool = True) -> bool:
    """디렉토리 존재 확인"""
    if path.exists() and path.is_dir():
        print(f"{Colors.GREEN}✓{Colors.NC} {name}: {path}")
        return True
    else:
        if critical:
            print(f"{Colors.RED}✗{Colors.NC} {name}: {path} (NOT FOUND)")
        else:
            print(f"{Colors.YELLOW}⚠{Colors.NC} {name}: {path} (OPTIONAL - NOT FOUND)")
        return False

def format_size(size_bytes: int) -> str:
    """바이트를 읽기 쉬운 형식으로 변환"""
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size_bytes < 1024.0:
            return f"{size_bytes:.1f}{unit}"
        size_bytes /= 1024.0
    return f"{size_bytes:.1f}PB"

def get_sample_dirs(aligned_dir: Path) -> List[str]:
    """aligned 디렉토리에서 샘플 목록 추출"""
    if not aligned_dir.exists():
        return []
    return [d.name for d in aligned_dir.iterdir() if d.is_dir()]

def verify_basic_structure() -> Tuple[int, int, int]:
    """기본 디렉토리 구조 검증"""
    print_section("1. 기본 디렉토리 구조 확인")
    
    results_dir = Path("results")
    pass_count = 0
    fail_count = 0
    warn_count = 0
    
    checks = [
        (results_dir, "Results 디렉토리", True),
        (results_dir / "trimmed", "Trimmed 디렉토리", True),
        (results_dir / "aligned", "Aligned 디렉토리", True),
        (results_dir / "counts", "Counts 디렉토리", True),
        (results_dir / "qc", "QC 디렉토리", False),
        (Path("logs"), "Logs 디렉토리", True),
    ]
    
    for path, name, critical in checks:
        if check_dir(path, name, critical):
            pass_count += 1
        elif critical:
            fail_count += 1
        else:
            warn_count += 1
    
    return pass_count, fail_count, warn_count

def verify_sample_files() -> Tuple[int, int, int]:
    """샘플별 파일 검증"""
    print_section("2. 샘플별 결과 파일 검증")
    
    pass_count = 0
    fail_count = 0
    warn_count = 0
    
    # Trimmed FASTQ 파일
    trimmed_dir = Path("results/trimmed")
    if trimmed_dir.exists():
        trimmed_files = list(trimmed_dir.glob("*.fastq.gz"))
        print(f"\n📁 Trimmed FASTQ 파일: {len(trimmed_files)}개")
        if trimmed_files:
            print(f"{Colors.GREEN}✓{Colors.NC} Trimmed FASTQ 파일 발견")
            pass_count += 1
            for f in trimmed_files[:5]:  # 처음 5개만 표시
                size = format_size(f.stat().st_size)
                print(f"  - {f.name} ({size})")
            if len(trimmed_files) > 5:
                print(f"  ... (총 {len(trimmed_files)}개)")
        else:
            print(f"{Colors.RED}✗{Colors.NC} Trimmed FASTQ 파일 없음")
            fail_count += 1
    
    # Aligned BAM 파일
    print("\n📁 Aligned BAM 파일:")
    aligned_dir = Path("results/aligned")
    samples = get_sample_dirs(aligned_dir)
    
    if samples:
        print(f"{Colors.GREEN}✓{Colors.NC} {len(samples)}개 샘플 발견")
        pass_count += 1
        
        for sample in samples:
            sample_dir = aligned_dir / sample
            bam_file = sample_dir / "Aligned.sortedByCoord.out.bam"
            log_file = sample_dir / "Log.final.out"
            
            if bam_file.exists():
                size = format_size(bam_file.stat().st_size)
                print(f"  {Colors.GREEN}✓{Colors.NC} {sample}: {bam_file.name} ({size})")
                pass_count += 1
            else:
                print(f"  {Colors.RED}✗{Colors.NC} {sample}: BAM 파일 없음")
                fail_count += 1
            
            if log_file.exists():
                print(f"    {Colors.GREEN}✓{Colors.NC} Log.final.out 존재")
                pass_count += 1
            else:
                print(f"    {Colors.YELLOW}⚠{Colors.NC} Log.final.out 없음")
                warn_count += 1
    else:
        print(f"{Colors.RED}✗{Colors.NC} BAM 파일 없음")
        fail_count += 1
    
    return pass_count, fail_count, warn_count

def verify_counts_matrix() -> Tuple[int, int, int]:
    """Counts matrix 검증"""
    print_section("3. Counts Matrix 검증")
    
    pass_count = 0
    fail_count = 0
    warn_count = 0
    
    counts_dir = Path("results/counts")
    raw_counts = counts_dir / "counts_matrix.txt"
    clean_counts = counts_dir / "counts_matrix_clean.csv"
    summary_file = counts_dir / "counts_matrix.txt.summary"
    
    # Raw counts
    if check_file(raw_counts, "Raw counts matrix", True):
        pass_count += 1
    else:
        fail_count += 1
    
    # Clean counts (CSV)
    if clean_counts.exists():
        if check_file(clean_counts, "Clean counts matrix (CSV)", True):
            pass_count += 1
            
            # CSV 내용 분석
            try:
                with open(clean_counts, 'r') as f:
                    lines = f.readlines()
                
                if lines:
                    # 헤더 (샘플 이름)
                    header = lines[0].strip().split(',')
                    sample_count = len(header) - 1  # gene_id 제외
                    gene_count = len(lines) - 1  # 헤더 제외
                    
                    print(f"\n  📊 Counts matrix 정보:")
                    print(f"  {Colors.GREEN}✓{Colors.NC} 유전자 수: {gene_count:,}")
                    print(f"  {Colors.GREEN}✓{Colors.NC} 샘플 수: {sample_count}")
                    pass_count += 2
                    
                    # 미리보기
                    print(f"\n  📊 처음 5줄:")
                    for line in lines[:5]:
                        print(f"    {line.strip()}")
            except Exception as e:
                print(f"{Colors.RED}✗{Colors.NC} CSV 파일 읽기 오류: {e}")
                fail_count += 1
        else:
            fail_count += 1
    
    # Summary 파일
    if summary_file.exists():
        if check_file(summary_file, "featureCounts summary", True):
            pass_count += 1
            
            print("\n  📊 featureCounts summary:")
            with open(summary_file, 'r') as f:
                print(f"    {f.read()}")
        else:
            fail_count += 1
    
    return pass_count, fail_count, warn_count

def verify_qc_reports() -> Tuple[int, int, int]:
    """QC 리포트 검증"""
    print_section("4. QC 리포트 검증")
    
    pass_count = 0
    fail_count = 0
    warn_count = 0
    
    # HTML QC 리포트
    qc_html = Path("results/qc_report.html")
    if qc_html.exists():
        if check_file(qc_html, "QC HTML 리포트", False):
            pass_count += 1
            
            # 샘플 개수 확인
            with open(qc_html, 'r', encoding='utf-8') as f:
                content = f.read()
                sample_count = content.count('sample-row')
                if sample_count > 0:
                    print(f"  {Colors.GREEN}✓{Colors.NC} 리포트 내 샘플 수: {sample_count}")
                    pass_count += 1
    else:
        print(f"{Colors.YELLOW}⚠{Colors.NC} QC HTML 리포트 없음 (선택사항)")
        warn_count += 1
    
    # QC JSON summary
    qc_json = Path("results/qc/qc_summary.json")
    if qc_json.exists():
        if check_file(qc_json, "QC JSON summary", False):
            pass_count += 1
            
            try:
                with open(qc_json, 'r') as f:
                    data = json.load(f)
                
                print("\n  📊 QC Summary 내용:")
                print(json.dumps(data, indent=2)[:500] + "...")
            except Exception as e:
                print(f"  {Colors.YELLOW}⚠{Colors.NC} JSON 파싱 오류: {e}")
                warn_count += 1
    else:
        print(f"{Colors.YELLOW}⚠{Colors.NC} QC JSON summary 없음 (선택사항)")
        warn_count += 1
    
    # FastQC 평가 결과
    fastqc_eval = Path("results/qc/fastqc_evaluation.json")
    if fastqc_eval.exists():
        if check_file(fastqc_eval, "FastQC 자동 평가 결과", False):
            pass_count += 1
            
            try:
                with open(fastqc_eval, 'r') as f:
                    data = json.load(f)
                
                print("\n  📊 FastQC 평가 요약:")
                if "overall_status" in data:
                    status = data["overall_status"]
                    emoji = "✓" if status == "PASS" else ("⚠" if status == "WARN" else "✗")
                    print(f"    {emoji} 전체 상태: {status}")
                
                if "samples" in data:
                    total = len(data["samples"])
                    passed = sum(1 for s in data["samples"].values() if s.get("status") == "PASS")
                    warned = sum(1 for s in data["samples"].values() if s.get("status") == "WARN")
                    failed = sum(1 for s in data["samples"].values() if s.get("status") == "FAIL")
                    
                    print(f"    총 샘플: {total}")
                    print(f"    PASS: {passed}, WARN: {warned}, FAIL: {failed}")
            except Exception as e:
                print(f"  {Colors.YELLOW}⚠{Colors.NC} JSON 파싱 오류: {e}")
                warn_count += 1
    
    return pass_count, fail_count, warn_count

def verify_logs() -> Tuple[int, int, int]:
    """로그 파일 검증"""
    print_section("5. 로그 파일 검증")
    
    pass_count = 0
    fail_count = 0
    warn_count = 0
    
    logs_dir = Path("logs")
    
    # 각 단계별 로그 개수
    cutadapt_logs = len(list((logs_dir / "cutadapt").glob("*.log"))) if (logs_dir / "cutadapt").exists() else 0
    star_logs = len(list((logs_dir / "star").glob("*.log"))) if (logs_dir / "star").exists() else 0
    fastqc_logs = len(list((logs_dir / "fastqc").glob("*.log"))) if (logs_dir / "fastqc").exists() else 0
    
    print(f"Cutadapt 로그: {cutadapt_logs}개")
    print(f"STAR 로그: {star_logs}개")
    print(f"FastQC 로그: {fastqc_logs}개")
    
    if cutadapt_logs > 0:
        pass_count += 1
    if star_logs > 0:
        pass_count += 1
    if fastqc_logs > 0:
        pass_count += 1
    
    # 주요 로그 파일
    if (logs_dir / "featurecounts.log").exists():
        print(f"{Colors.GREEN}✓{Colors.NC} featureCounts 로그 존재")
        pass_count += 1
    
    if (logs_dir / "qc_report.log").exists():
        print(f"{Colors.GREEN}✓{Colors.NC} QC 리포트 로그 존재")
        pass_count += 1
    
    return pass_count, fail_count, warn_count

def show_disk_usage():
    """디스크 사용량 표시"""
    print_section("6. 디스크 사용량")
    
    results_dir = Path("results")
    if not results_dir.exists():
        print(f"{Colors.RED}✗{Colors.NC} results 디렉토리 없음")
        return
    
    # 디렉토리별 용량
    print("📊 디렉토리별 용량:")
    subdirs = [d for d in results_dir.iterdir() if d.is_dir()]
    
    sizes = []
    for subdir in subdirs:
        total_size = sum(f.stat().st_size for f in subdir.rglob('*') if f.is_file())
        sizes.append((subdir.name, total_size))
    
    sizes.sort(key=lambda x: x[1], reverse=True)
    
    for name, size in sizes:
        print(f"  {name}: {format_size(size)}")
    
    # 총 용량
    total_size = sum(f.stat().st_size for f in results_dir.rglob('*') if f.is_file())
    print(f"\n📊 총 용량: {format_size(total_size)}")

def main():
    """메인 함수"""
    print_header("RNA-seq Pipeline 결과 구조 검증")
    
    total_pass = 0
    total_fail = 0
    total_warn = 0
    
    # 각 검증 단계 실행
    p, f, w = verify_basic_structure()
    total_pass += p
    total_fail += f
    total_warn += w
    
    p, f, w = verify_sample_files()
    total_pass += p
    total_fail += f
    total_warn += w
    
    p, f, w = verify_counts_matrix()
    total_pass += p
    total_fail += f
    total_warn += w
    
    p, f, w = verify_qc_reports()
    total_pass += p
    total_fail += f
    total_warn += w
    
    p, f, w = verify_logs()
    total_pass += p
    total_fail += f
    total_warn += w
    
    # 디스크 사용량
    show_disk_usage()
    
    # 최종 요약
    print_header("검증 결과 요약")
    print(f"{Colors.GREEN}✓ PASS:{Colors.NC} {total_pass}")
    print(f"{Colors.YELLOW}⚠ WARN:{Colors.NC} {total_warn}")
    print(f"{Colors.RED}✗ FAIL:{Colors.NC} {total_fail}")
    
    total = total_pass + total_warn + total_fail
    if total > 0:
        success_rate = (total_pass * 100) // total
        print(f"\n성공률: {success_rate}%")
    
    print()
    
    if total_fail == 0:
        print(f"{Colors.GREEN}🎉 모든 필수 검증 통과!{Colors.NC}")
        if total_warn > 0:
            print(f"{Colors.YELLOW}⚠️  {total_warn}개의 선택사항 파일이 누락되었습니다.{Colors.NC}")
        print()
        print("✅ 브랜치 병합(merge) 준비 완료")
        return 0
    else:
        print(f"{Colors.RED}❌ {total_fail}개의 필수 검증 실패{Colors.NC}")
        print()
        print("⚠️  위 문제를 해결한 후 병합하세요.")
        return 1

if __name__ == "__main__":
    sys.exit(main())
