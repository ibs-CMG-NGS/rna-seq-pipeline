#!/usr/bin/env python3
"""
Sample Sheet Converter

마스터 samples.csv를 각 파이프라인 형식으로 변환하는 스크립트
"""

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import List, Dict


class SampleSheetConverter:
    """샘플 시트 변환기"""
    
    # 마스터 시트 필수 컬럼
    REQUIRED_COLUMNS = [
        'project_id', 'sample_id', 'sample_name', 'condition', 
        'replicate', 'library_type', 'read_type', 'fastq_1', 
        'species', 'genome_build'
    ]
    
    def __init__(self, master_csv_path: str):
        """
        Args:
            master_csv_path: 마스터 samples.csv 파일 경로
        """
        self.master_path = Path(master_csv_path)
        self.samples = self._load_master_sheet()
        
    def _load_master_sheet(self) -> List[Dict]:
        """마스터 샘플 시트 로드 및 검증"""
        if not self.master_path.exists():
            raise FileNotFoundError(f"Master sample sheet not found: {self.master_path}")
        
        samples = []
        with open(self.master_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            
            # 필수 컬럼 확인
            missing_cols = set(self.REQUIRED_COLUMNS) - set(reader.fieldnames)
            if missing_cols:
                raise ValueError(f"Missing required columns: {missing_cols}")
            
            for row in reader:
                # 빈 행 건너뛰기
                if not row.get('sample_id'):
                    continue
                    
                # FASTQ 파일 존재 확인 (선택사항)
                # if not Path(row['fastq_1']).exists():
                #     print(f"Warning: FASTQ not found: {row['fastq_1']}", file=sys.stderr)
                
                samples.append(row)
        
        return samples
    
    def filter_by_library(self, library_type: str) -> List[Dict]:
        """라이브러리 타입으로 샘플 필터링"""
        return [s for s in self.samples if s['library_type'].lower() == library_type.lower()]
    
    def filter_by_project(self, project_id: str) -> List[Dict]:
        """프로젝트 ID로 샘플 필터링"""
        return [s for s in self.samples if s['project_id'] == project_id]
    
    def to_snakemake_tsv(self, output_path: str, library_type: str = 'RNA-seq'):
        """
        Snakemake용 TSV 형식으로 변환
        
        Args:
            output_path: 출력 파일 경로
            library_type: 라이브러리 타입 (기본: RNA-seq)
        """
        samples = self.filter_by_library(library_type)
        
        if not samples:
            raise ValueError(f"No samples found for library type: {library_type}")
        
        with open(output_path, 'w', encoding='utf-8', newline='') as f:
            # TSV 헤더
            f.write("sample_id\tcondition\treplicate\tfastq_r1\tfastq_r2\tnotes\n")
            
            for sample in samples:
                row = [
                    sample['sample_id'],
                    sample['condition'],
                    sample['replicate'],
                    sample['fastq_1'],
                    sample.get('fastq_2', ''),  # SE인 경우 빈 값
                    sample.get('notes', '')
                ]
                f.write('\t'.join(row) + '\n')
        
        print(f"✓ Snakemake TSV generated: {output_path} ({len(samples)} samples)")
    
    def to_wdl_json(self, output_path: str, library_type: str = 'WGS'):
        """
        WDL용 JSON 형식으로 변환
        
        Args:
            output_path: 출력 파일 경로
            library_type: 라이브러리 타입 (기본: WGS)
        """
        samples = self.filter_by_library(library_type)
        
        if not samples:
            raise ValueError(f"No samples found for library type: {library_type}")
        
        wdl_input = {
            "workflow.samples": []
        }
        
        for sample in samples:
            sample_dict = {
                "sample_id": sample['sample_id'],
                "fastq_1": sample['fastq_1'],
                "read_type": sample['read_type']
            }
            
            if sample.get('fastq_2'):
                sample_dict["fastq_2"] = sample['fastq_2']
            
            if sample.get('genome_build'):
                sample_dict["reference_genome"] = sample['genome_build']
            
            wdl_input["workflow.samples"].append(sample_dict)
        
        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(wdl_input, f, indent=2)
        
        print(f"✓ WDL JSON generated: {output_path} ({len(samples)} samples)")
    
    def to_nextflow_csv(self, output_path: str, library_type: str = 'ATAC-seq'):
        """
        Nextflow용 CSV 형식으로 변환
        
        Args:
            output_path: 출력 파일 경로
            library_type: 라이브러리 타입 (기본: ATAC-seq)
        """
        samples = self.filter_by_library(library_type)
        
        if not samples:
            raise ValueError(f"No samples found for library type: {library_type}")
        
        with open(output_path, 'w', encoding='utf-8', newline='') as f:
            writer = csv.writer(f)
            
            # CSV 헤더
            writer.writerow(['sample', 'fastq_1', 'fastq_2', 'replicate'])
            
            for sample in samples:
                row = [
                    sample['sample_id'],
                    sample['fastq_1'],
                    sample.get('fastq_2', ''),
                    sample['replicate']
                ]
                writer.writerow(row)
        
        print(f"✓ Nextflow CSV generated: {output_path} ({len(samples)} samples)")
    
    def generate_all(self, output_dir: str, project_id: str = None):
        """
        모든 파이프라인 형식으로 변환
        
        Args:
            output_dir: 출력 디렉토리
            project_id: 프로젝트 ID (선택사항, 지정시 해당 프로젝트만)
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        samples = self.samples
        if project_id:
            samples = self.filter_by_project(project_id)
        
        # 라이브러리 타입별로 샘플 그룹화
        library_types = set(s['library_type'] for s in samples)
        
        print(f"\n{'='*70}")
        print(f"Converting master sample sheet: {self.master_path}")
        print(f"Output directory: {output_path}")
        print(f"Total samples: {len(samples)}")
        print(f"Library types: {', '.join(library_types)}")
        print(f"{'='*70}\n")
        
        # RNA-seq → Snakemake TSV
        if any('rna' in lt.lower() for lt in library_types):
            try:
                self.to_snakemake_tsv(
                    str(output_path / 'samples_rnaseq.tsv'),
                    'RNA-seq'
                )
            except ValueError as e:
                print(f"⚠ Skipping RNA-seq: {e}")
        
        # WGS → WDL JSON
        if 'WGS' in library_types:
            try:
                self.to_wdl_json(
                    str(output_path / 'inputs_wgs.json'),
                    'WGS'
                )
            except ValueError as e:
                print(f"⚠ Skipping WGS: {e}")
        
        # ATAC-seq → Nextflow CSV
        if any('atac' in lt.lower() for lt in library_types):
            try:
                self.to_nextflow_csv(
                    str(output_path / 'samplesheet_atacseq.csv'),
                    'ATAC-seq'
                )
            except ValueError as e:
                print(f"⚠ Skipping ATAC-seq: {e}")
        
        print(f"\n✅ Conversion complete!\n")


def create_template(output_path: str):
    """마스터 샘플 시트 템플릿 생성"""
    template_data = [
        {
            'project_id': 'PROJECT_2025',
            'sample_id': 'Sample_Ctrl_1',
            'sample_name': 'Control_Rep1',
            'condition': 'Control',
            'replicate': '1',
            'sequencing_platform': 'Illumina',
            'library_type': 'RNA-seq',
            'read_type': 'paired-end',
            'fastq_1': '/path/to/Ctrl_1_R1.fastq.gz',
            'fastq_2': '/path/to/Ctrl_1_R2.fastq.gz',
            'species': 'human',
            'genome_build': 'GRCh38',
            'batch': '1',
            'sex': 'F',
            'tissue': 'fibroblast',
            'treatment': 'none',
            'time_point': '0h',
            'notes': 'Control sample replicate 1'
        },
        {
            'project_id': 'PROJECT_2025',
            'sample_id': 'Sample_Treat_1',
            'sample_name': 'Treatment_Rep1',
            'condition': 'Treatment',
            'replicate': '1',
            'sequencing_platform': 'Illumina',
            'library_type': 'RNA-seq',
            'read_type': 'paired-end',
            'fastq_1': '/path/to/Treat_1_R1.fastq.gz',
            'fastq_2': '/path/to/Treat_1_R2.fastq.gz',
            'species': 'human',
            'genome_build': 'GRCh38',
            'batch': '1',
            'sex': 'F',
            'tissue': 'fibroblast',
            'treatment': 'drug_100uM',
            'time_point': '24h',
            'notes': 'Treatment sample replicate 1'
        }
    ]
    
    with open(output_path, 'w', encoding='utf-8', newline='') as f:
        if template_data:
            writer = csv.DictWriter(f, fieldnames=template_data[0].keys())
            writer.writeheader()
            writer.writerows(template_data)
    
    print(f"✓ Template created: {output_path}")
    print(f"  Edit this file and fill in your sample information")


def main():
    parser = argparse.ArgumentParser(
        description='Convert master sample sheet to pipeline-specific formats',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Create a template
  python convert_sample_sheet.py --create-template samples.csv
  
  # Convert for all pipelines
  python convert_sample_sheet.py samples.csv -o config/
  
  # Convert only for specific project
  python convert_sample_sheet.py samples.csv -o config/ --project H2O2_2025
  
  # Convert to Snakemake TSV only
  python convert_sample_sheet.py samples.csv --snakemake samples.tsv
        """
    )
    
    parser.add_argument(
        'input',
        nargs='?',
        help='Master samples.csv file'
    )
    parser.add_argument(
        '-o', '--output-dir',
        help='Output directory for all pipeline formats'
    )
    parser.add_argument(
        '--project',
        help='Filter by project ID'
    )
    parser.add_argument(
        '--snakemake',
        metavar='FILE',
        help='Generate Snakemake TSV format'
    )
    parser.add_argument(
        '--wdl',
        metavar='FILE',
        help='Generate WDL JSON format'
    )
    parser.add_argument(
        '--nextflow',
        metavar='FILE',
        help='Generate Nextflow CSV format'
    )
    parser.add_argument(
        '--create-template',
        metavar='FILE',
        help='Create a template samples.csv file'
    )
    
    args = parser.parse_args()
    
    # 템플릿 생성
    if args.create_template:
        create_template(args.create_template)
        return 0
    
    # 입력 파일 확인
    if not args.input:
        parser.error("Input file required (or use --create-template)")
    
    try:
        converter = SampleSheetConverter(args.input)
        
        # 모든 형식 생성
        if args.output_dir:
            converter.generate_all(args.output_dir, args.project)
        
        # 개별 형식 생성
        if args.snakemake:
            converter.to_snakemake_tsv(args.snakemake)
        
        if args.wdl:
            converter.to_wdl_json(args.wdl)
        
        if args.nextflow:
            converter.to_nextflow_csv(args.nextflow)
        
        if not any([args.output_dir, args.snakemake, args.wdl, args.nextflow]):
            parser.error("Specify at least one output option (-o, --snakemake, --wdl, or --nextflow)")
        
        return 0
        
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
