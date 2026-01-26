#!/usr/bin/env python3
"""
FastQC 자동 평가 스크립트 테스트

이 스크립트는 evaluate_fastqc.py가 올바르게 작동하는지 간단히 테스트합니다.
"""

import os
import sys
import tempfile
import shutil
from pathlib import Path

# 테스트용 fastqc_data.txt 생성
SAMPLE_FASTQC_DATA_GOOD = """##FastQC	0.12.1
>>Basic Statistics	pass
#Measure	Value
Filename	Sample_Good_1.fastq.gz
File type	Conventional base calls
Encoding	Sanger / Illumina 1.9
Total Sequences	5280544
Sequences flagged as poor quality	0
Sequence length	101
%GC	48
>>END_MODULE
>>Per base sequence quality	pass
#Base	Mean	Median	Lower Quartile	Upper Quartile	10th Percentile	90th Percentile
1	37.5	38.0	36.0	39.0	35.0	39.0
50	35.2	36.0	34.0	37.0	33.0	37.0
101	34.8	35.0	33.0	36.0	32.0	36.0
>>END_MODULE
>>Per tile sequence quality	warn
>>END_MODULE
>>Per sequence quality scores	pass
#Quality	Count
30	1000000.0
35	3000000.0
38	1280544.0
>>END_MODULE
>>Per base sequence content	fail
>>END_MODULE
>>Per sequence GC content	pass
>>END_MODULE
>>Per base N content	pass
#Base	N-Count
1	0.1
50	0.2
101	0.3
>>END_MODULE
>>Sequence Length Distribution	pass
>>END_MODULE
>>Sequence Duplication Levels	fail
>>END_MODULE
>>Overrepresented sequences	fail
>>END_MODULE
>>Adapter Content	warn
#Position	Illumina Universal Adapter
1	0.1
50	2.5
101	5.8
>>END_MODULE
"""

SAMPLE_FASTQC_DATA_BAD = """##FastQC	0.12.1
>>Basic Statistics	fail
#Measure	Value
Filename	Sample_Bad_1.fastq.gz
File type	Conventional base calls
Encoding	Sanger / Illumina 1.9
Total Sequences	856432
Sequences flagged as poor quality	0
Sequence length	101
%GC	28
>>END_MODULE
>>Per base sequence quality	fail
#Base	Mean	Median	Lower Quartile	Upper Quartile	10th Percentile	90th Percentile
1	25.5	26.0	22.0	28.0	20.0	29.0
50	18.2	19.0	15.0	21.0	13.0	22.0
101	15.8	16.0	12.0	18.0	10.0	19.0
>>END_MODULE
>>Per tile sequence quality	warn
>>END_MODULE
>>Per sequence quality scores	fail
#Quality	Count
15	300000.0
20	400000.0
25	156432.0
>>END_MODULE
>>Per base sequence content	fail
>>END_MODULE
>>Per sequence GC content	pass
>>END_MODULE
>>Per base N content	fail
#Base	N-Count
1	2.1
50	8.5
101	12.3
>>END_MODULE
>>Sequence Length Distribution	pass
>>END_MODULE
>>Sequence Duplication Levels	fail
>>END_MODULE
>>Overrepresented sequences	fail
>>END_MODULE
>>Adapter Content	fail
#Position	Illumina Universal Adapter
1	1.5
50	15.8
101	28.3
>>END_MODULE
"""

def create_mock_fastqc_result(output_dir: Path, sample_name: str, data_content: str):
    """모의 FastQC 결과 생성"""
    # ZIP 파일 (실제로는 생성 안하고 디렉토리만)
    zip_path = output_dir / f"{sample_name}_fastqc.zip"
    zip_path.touch()
    
    # 압축 해제된 디렉토리
    data_dir = output_dir / f"{sample_name}_fastqc"
    data_dir.mkdir(exist_ok=True)
    
    # fastqc_data.txt 생성
    data_file = data_dir / "fastqc_data.txt"
    with open(data_file, 'w') as f:
        f.write(data_content)
    
    print(f"✓ Created mock FastQC result: {sample_name}")

def run_test():
    """테스트 실행"""
    print("=" * 80)
    print("FastQC Auto-Evaluation Test")
    print("=" * 80)
    print()
    
    # 임시 디렉토리 생성
    with tempfile.TemporaryDirectory() as tmpdir:
        qc_dir = Path(tmpdir) / "qc"
        qc_dir.mkdir()
        
        print("1. Creating mock FastQC results...")
        print()
        
        # 좋은 샘플 3개
        for i in range(1, 4):
            create_mock_fastqc_result(
                qc_dir, 
                f"Sample_Good_{i}", 
                SAMPLE_FASTQC_DATA_GOOD
            )
        
        # 나쁜 샘플 2개
        for i in range(1, 3):
            create_mock_fastqc_result(
                qc_dir, 
                f"Sample_Bad_{i}", 
                SAMPLE_FASTQC_DATA_BAD
            )
        
        print()
        print("2. Running evaluation...")
        print()
        
        # 평가 실행
        report_path = qc_dir / "evaluation_report.txt"
        json_path = qc_dir / "evaluation_report.json"
        
        cmd = [
            sys.executable,
            "src/evaluate_fastqc.py",
            str(qc_dir),
            "-o", str(report_path),
            "--json", str(json_path)
        ]
        
        import subprocess
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        print(result.stdout)
        
        if result.returncode != 0:
            print("ERROR:", result.stderr)
            return False
        
        print()
        print("3. Checking results...")
        print()
        
        # 리포트 파일 확인
        if report_path.exists():
            print(f"✓ Report generated: {report_path}")
            print()
            print("=" * 80)
            print("Report Preview:")
            print("=" * 80)
            with open(report_path, 'r', encoding='utf-8') as f:
                lines = f.readlines()
                # 처음 50줄만 출력
                for line in lines[:50]:
                    print(line.rstrip())
                if len(lines) > 50:
                    print(f"\n... ({len(lines) - 50} more lines)")
        else:
            print("✗ Report not generated")
            return False
        
        if json_path.exists():
            print()
            print(f"✓ JSON generated: {json_path}")
            
            # JSON 내용 확인
            import json
            with open(json_path, 'r') as f:
                data = json.load(f)
                print(f"  - Total samples: {len(data)}")
                status_count = {}
                for item in data:
                    status = item['status']
                    status_count[status] = status_count.get(status, 0) + 1
                for status, count in sorted(status_count.items()):
                    print(f"  - {status}: {count}")
        else:
            print("✗ JSON not generated")
            return False
        
        print()
        print("=" * 80)
        print("✅ Test completed successfully!")
        print("=" * 80)
        
        return True

if __name__ == '__main__':
    success = run_test()
    sys.exit(0 if success else 1)
