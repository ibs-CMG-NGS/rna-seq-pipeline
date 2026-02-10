#!/usr/bin/env python3
"""
Manifest Generator

파이프라인 실행 결과에서 manifest.json 파일을 생성합니다.
manifest.json은 에이전트가 최종 결과물을 식별하고 검증하는데 사용됩니다.
"""

import argparse
import json
import hashlib
import os
import sys
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False


class ManifestGenerator:
    """Manifest 파일 생성기"""
    
    def __init__(self, sample_dir: str, sample_id: str, project_id: str, pipeline_type: str = "rna-seq", 
                 sample_metadata: Optional[Dict] = None):
        """
        Args:
            sample_dir: 샘플 디렉토리 경로
            sample_id: 샘플 ID
            project_id: 프로젝트 ID
            pipeline_type: 파이프라인 타입
            sample_metadata: 샘플 메타데이터 (condition, replicate 등)
        """
        self.sample_dir = Path(sample_dir)
        self.sample_id = sample_id
        self.project_id = project_id
        self.pipeline_type = pipeline_type
        self.sample_metadata = sample_metadata or {}
        
        self.final_outputs_dir = self.sample_dir / "final_outputs"
        self.intermediate_dir = self.sample_dir / "intermediate"
        self.metadata_dir = self.sample_dir / "metadata"
        
    def calculate_md5(self, file_path: Path, chunk_size: int = 8192) -> str:
        """파일의 MD5 체크섬 계산"""
        md5 = hashlib.md5()
        
        try:
            with open(file_path, 'rb') as f:
                while chunk := f.read(chunk_size):
                    md5.update(chunk)
            return md5.hexdigest()
        except Exception as e:
            print(f"Warning: Could not calculate MD5 for {file_path}: {e}", file=sys.stderr)
            return "N/A"
    
    def get_file_info(self, file_path: Path, description: str = "") -> Dict:
        """파일 정보를 딕셔너리로 반환"""
        if not file_path.exists():
            return None
        
        # final_outputs 디렉토리로부터의 상대 경로
        try:
            rel_path = file_path.relative_to(self.final_outputs_dir)
        except ValueError:
            rel_path = file_path.name
        
        return {
            "path": str(rel_path),
            "absolute_path": str(file_path.absolute()),
            "md5": self.calculate_md5(file_path),
            "size_bytes": file_path.stat().st_size,
            "modified": datetime.fromtimestamp(file_path.stat().st_mtime).isoformat(),
            "description": description
        }
    
    def find_rna_seq_outputs(self) -> Dict:
        """RNA-seq 파이프라인의 최종 결과물 찾기"""
        outputs = {}
        
        # BAM 파일
        bam_dir = self.final_outputs_dir / "bam"
        if bam_dir.exists():
            bam_file = bam_dir / "aligned.sorted.bam"
            bai_file = bam_dir / "aligned.sorted.bam.bai"
            
            if bam_file.exists():
                outputs["aligned_bam"] = self.get_file_info(
                    bam_file, 
                    "STAR aligned reads, sorted by coordinate"
                )
            
            if bai_file.exists():
                outputs["bam_index"] = self.get_file_info(
                    bai_file,
                    "BAM index file for IGV/samtools"
                )
        
        # Gene counts
        counts_dir = self.final_outputs_dir / "counts"
        if counts_dir.exists():
            counts_file = counts_dir / "gene_counts.csv"
            if counts_file.exists():
                outputs["gene_counts"] = self.get_file_info(
                    counts_file,
                    "Gene-level read counts from featureCounts"
                )
        
        # QC 리포트
        qc_dir = self.final_outputs_dir / "qc"
        if qc_dir.exists():
            multiqc_file = qc_dir / "multiqc_report.html"
            qc_summary_file = qc_dir / "qc_summary.json"
            
            if multiqc_file.exists():
                outputs["qc_report"] = self.get_file_info(
                    multiqc_file,
                    "Comprehensive QC report from MultiQC"
                )
            
            if qc_summary_file.exists():
                outputs["qc_summary"] = self.get_file_info(
                    qc_summary_file,
                    "Machine-readable QC metrics"
                )
        
        return outputs
    
    def extract_qc_metrics(self) -> Dict:
        """QC summary에서 주요 메트릭 추출"""
        qc_summary_file = self.final_outputs_dir / "qc" / "qc_summary.json"
        
        if qc_summary_file.exists():
            try:
                with open(qc_summary_file, 'r') as f:
                    return json.load(f)
            except Exception as e:
                print(f"Warning: Could not read QC summary: {e}", file=sys.stderr)
        
        return {
            "overall_status": "UNKNOWN",
            "note": "QC summary not found or could not be parsed"
        }
    
    def suggest_next_steps(self, qc_metrics: Dict) -> List[str]:
        """QC 상태에 따라 다음 단계 제안"""
        status = qc_metrics.get("overall_status", "UNKNOWN")
        
        if status == "PASS":
            return [
                "differential_expression",
                "pathway_analysis",
                "gene_set_enrichment"
            ]
        elif status == "WARN":
            return [
                "review_qc_report",
                "conditional_differential_expression"
            ]
        elif status == "FAIL":
            return [
                "exclude_from_analysis",
                "consider_resequencing"
            ]
        else:
            return ["review_outputs"]
    
    def _extract_metadata(self) -> Dict:
        """샘플 메타데이터 추출 및 정리"""
        if not self.sample_metadata:
            return {}
        
        # 관심 있는 메타데이터 필드만 추출
        relevant_fields = [
            'sample_name', 'condition', 'replicate', 'sequencing_platform',
            'library_type', 'read_type', 'species', 'genome_build',
            'batch', 'tissue', 'treatment', 'time_point', 'notes'
        ]
        
        metadata = {}
        for field in relevant_fields:
            if field in self.sample_metadata and self.sample_metadata[field]:
                value = self.sample_metadata[field]
                # pandas가 있으면 NaN 체크, 없으면 None과 빈 문자열만 체크
                if HAS_PANDAS:
                    if pd.notna(value) and value != '':
                        metadata[field] = value
                else:
                    if value is not None and value != '':
                        metadata[field] = value
        
        return metadata
    
    def generate_manifest(self, output_path: Optional[str] = None) -> Dict:
        """Manifest 생성"""
        # 최종 결과물 찾기
        if self.pipeline_type == "rna-seq":
            final_outputs = self.find_rna_seq_outputs()
        else:
            raise ValueError(f"Unsupported pipeline type: {self.pipeline_type}")
        
        # QC 메트릭 추출
        qc_metrics = self.extract_qc_metrics()
        
        # Manifest 구성
        manifest = {
            "sample_id": self.sample_id,
            "project_id": self.project_id,
            "pipeline_type": self.pipeline_type,
            "pipeline_version": "1.0.0",  # TODO: 자동으로 가져오기
            "execution_date": datetime.now().isoformat(),
            "status": "completed" if final_outputs else "incomplete",
            
            # 샘플 메타데이터 추가
            "sample_metadata": self._extract_metadata(),
            
            "final_outputs": final_outputs,
            
            "qc_metrics": qc_metrics,
            
            "next_steps": self.suggest_next_steps(qc_metrics),
            
            "directory_structure": {
                "sample_dir": str(self.sample_dir.absolute()),
                "final_outputs": str(self.final_outputs_dir.absolute()),
                "intermediate": str(self.intermediate_dir.absolute()),
                "metadata": str(self.metadata_dir.absolute())
            }
        }
        
        # 파일로 저장
        if output_path is None:
            output_path = self.final_outputs_dir / "manifest.json"
        
        output_file = Path(output_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(manifest, f, indent=2, ensure_ascii=False)
        
        return manifest
    
    def validate_manifest(self, manifest_path: str) -> bool:
        """Manifest 파일 검증 (MD5 체크섬)"""
        with open(manifest_path, 'r') as f:
            manifest = json.load(f)
        
        all_valid = True
        
        for output_name, output_info in manifest.get("final_outputs", {}).items():
            file_path = Path(output_info["absolute_path"])
            
            if not file_path.exists():
                print(f"❌ MISSING: {output_name} - {file_path}", file=sys.stderr)
                all_valid = False
                continue
            
            expected_md5 = output_info.get("md5")
            if expected_md5 == "N/A":
                print(f"⚠️  SKIP: {output_name} - No MD5 available")
                continue
            
            actual_md5 = self.calculate_md5(file_path)
            
            if actual_md5 == expected_md5:
                print(f"✅ OK: {output_name}")
            else:
                print(f"❌ CHECKSUM MISMATCH: {output_name}", file=sys.stderr)
                print(f"   Expected: {expected_md5}", file=sys.stderr)
                print(f"   Actual: {actual_md5}", file=sys.stderr)
                all_valid = False
        
        return all_valid


def main():
    parser = argparse.ArgumentParser(
        description='Generate manifest.json for pipeline outputs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate manifest for a sample
  python generate_manifest.py \\
      --sample-dir /data/results/PROJECT/SAMPLE/rna-seq \\
      --sample-id SAMPLE \\
      --project-id PROJECT
  
  # Validate existing manifest
  python generate_manifest.py \\
      --validate /data/results/PROJECT/SAMPLE/rna-seq/final_outputs/manifest.json
        """
    )
    
    parser.add_argument(
        '--sample-dir',
        help='Sample directory path (e.g., /results/PROJECT/SAMPLE/rna-seq)'
    )
    parser.add_argument(
        '--sample-id',
        help='Sample ID'
    )
    parser.add_argument(
        '--project-id',
        help='Project ID'
    )
    parser.add_argument(
        '--pipeline-type',
        default='rna-seq',
        help='Pipeline type (default: rna-seq)'
    )
    parser.add_argument(
        '--sample-metadata',
        help='Sample metadata as JSON string (e.g., \'{"condition": "Control", "replicate": 1}\')'
    )
    parser.add_argument(
        '-o', '--output',
        help='Output manifest.json path (default: {sample_dir}/final_outputs/manifest.json)'
    )
    parser.add_argument(
        '--validate',
        metavar='MANIFEST',
        help='Validate an existing manifest file'
    )
    
    args = parser.parse_args()
    
    try:
        if args.validate:
            # Manifest 검증 모드
            manifest_path = Path(args.validate)
            sample_dir = manifest_path.parent.parent
            
            with open(manifest_path, 'r') as f:
                manifest = json.load(f)
            
            generator = ManifestGenerator(
                str(sample_dir),
                manifest['sample_id'],
                manifest['project_id'],
                manifest['pipeline_type']
            )
            
            print(f"\n{'='*70}")
            print(f"Validating Manifest: {manifest_path}")
            print(f"{'='*70}\n")
            
            is_valid = generator.validate_manifest(str(manifest_path))
            
            print(f"\n{'='*70}")
            if is_valid:
                print("✅ Manifest validation PASSED")
            else:
                print("❌ Manifest validation FAILED")
            print(f"{'='*70}\n")
            
            return 0 if is_valid else 1
            
        else:
            # Manifest 생성 모드
            if not all([args.sample_dir, args.sample_id, args.project_id]):
                parser.error("--sample-dir, --sample-id, and --project-id are required for generation")
            
            # Parse sample metadata if provided
            sample_metadata = None
            if args.sample_metadata:
                try:
                    sample_metadata = json.loads(args.sample_metadata)
                except json.JSONDecodeError as e:
                    parser.error(f"Invalid JSON in --sample-metadata: {e}")
            
            generator = ManifestGenerator(
                args.sample_dir,
                args.sample_id,
                args.project_id,
                args.pipeline_type,
                sample_metadata
            )
            
            print(f"\n{'='*70}")
            print(f"Generating Manifest")
            print(f"  Sample: {args.sample_id}")
            print(f"  Project: {args.project_id}")
            print(f"  Pipeline: {args.pipeline_type}")
            print(f"  Directory: {args.sample_dir}")
            print(f"{'='*70}\n")
            
            manifest = generator.generate_manifest(args.output)
            
            output_file = args.output or f"{args.sample_dir}/final_outputs/manifest.json"
            
            print(f"✅ Manifest generated: {output_file}")
            print(f"\nSummary:")
            print(f"  Status: {manifest['status']}")
            print(f"  Final outputs: {len(manifest['final_outputs'])}")
            print(f"  QC status: {manifest['qc_metrics'].get('overall_status', 'N/A')}")
            print(f"  Next steps: {', '.join(manifest['next_steps'])}\n")
            
            return 0
            
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
