# Batch Runner

여러 RNA-seq 프로젝트를 LLM 없이 순차적으로 무인 실행하는 스크립트.

## 파일 위치

```
scripts/batch_runner.py          # 실행 스크립트
config/batch/example_batch.yaml  # 배치 설정 템플릿
batch_logs/                      # 실행 로그 및 체크포인트 (자동 생성)
```

## 실행 흐름

프로젝트당 아래 순서로 실행되며, 실패 시 해당 프로젝트만 건너뛰고 다음 프로젝트로 진행한다.

```
1. validate_input_data   → FASTQ, 게놈, 디스크 공간 사전 확인
2. run_pipeline          → Snakemake 실행 (blocking, 멀티코어)
3. read_qc_results       → FastQC / MultiQC / STAR 매핑률 확인
4. run_bridge            → DE/GO 파이프라인 연결 (run_bridge: true 인 경우만)
```

## 사용법

```bash
# conda 환경 활성화 필수
conda activate rnaseq  # 또는 snakemake가 포함된 환경

# dry-run: 설정/경로 확인만 (Snakemake 실제 실행 없음)
python scripts/batch_runner.py config/batch/my_batch.yaml --dry-run

# 실제 실행
python scripts/batch_runner.py config/batch/my_batch.yaml --cores 12

# 중단된 실행 재개 (완료된 step은 자동 스킵)
python scripts/batch_runner.py config/batch/my_batch.yaml --cores 12 --resume
```

### 옵션

| 옵션 | 기본값 | 설명 |
|------|--------|------|
| `--dry-run` | false | Snakemake를 dry-run 모드로 실행 |
| `--resume` | false | `.done` 체크포인트가 있는 step 스킵 |
| `--cores N` | batch YAML의 `cores` 값 | Snakemake 코어 수 재정의 |

## 배치 설정 파일

```yaml
# config/batch/my_batch.yaml

pipeline_root: /data_3tb/shared/rna-seq-pipeline

# run_bridge: true 인 프로젝트가 있을 때만 필요
de_pipeline_dir: /home/user/ngs_pipeline/jinhyeong/RNA-Seq_DE_GO_analysis

cores: 12  # Snakemake 기본 코어 수

projects:
  - project_id: mouse-chd8
    config_file: config/projects/config_mouse_chd8_local.yaml
    run_bridge: true

  - project_id: human-project
    config_file: config/projects/config_human_local.yaml
    run_bridge: false
```

### 필드 설명

| 필드 | 필수 | 설명 |
|------|------|------|
| `pipeline_root` | 권장 | 파이프라인 루트 경로 (현재는 참고용) |
| `de_pipeline_dir` | 조건부 | `run_bridge: true`가 있으면 필수 |
| `cores` | 아니오 | 기본값 8, `--cores`로 재정의 가능 |
| `projects[].project_id` | 예 | 로그 디렉토리명 및 요약에 사용 |
| `projects[].config_file` | 예 | 파이프라인 config YAML 경로 |
| `projects[].run_bridge` | 아니오 | 기본값 false |

## 출력

### 체크포인트 파일

`batch_logs/<project_id>/<step>.done` 형식으로 완료된 step 추적.

```
batch_logs/
├── batch.log                   # 전체 실행 로그
├── summary.csv                 # 전체 결과 요약
├── mouse-chd8/
│   ├── runner.log              # 프로젝트별 상세 로그
│   ├── validate.done
│   ├── pipeline.done
│   ├── qc_check.done
│   └── bridge.done
└── human-project/
    ├── runner.log
    ├── validate.done
    └── pipeline.done           # bridge는 run_bridge: false이므로 없음
```

### summary.csv

```
project_id,status,failed_step,duration_sec,qc_pass,mapping_pct,error
mouse-chd8,success,,5432.1,True,92.4%,
human-project,failed,pipeline,12.3,,,Pipeline execution failed: ...
```

## 재실행 (resume)

중간에 실패하거나 중단된 경우 `--resume`으로 재실행하면 이미 완료된 step을 자동으로 건너뛴다.

```bash
# 실패한 프로젝트의 특정 step부터 재시작
python scripts/batch_runner.py config/batch/my_batch.yaml --cores 12 --resume
```

처음부터 다시 실행하려면 해당 프로젝트의 `.done` 파일을 삭제한다.

```bash
# 특정 프로젝트의 pipeline step부터 다시
rm batch_logs/mouse-chd8/pipeline.done batch_logs/mouse-chd8/qc_check.done batch_logs/mouse-chd8/bridge.done
python scripts/batch_runner.py config/batch/my_batch.yaml --cores 12 --resume

# 특정 프로젝트 전체 초기화
rm -rf batch_logs/mouse-chd8/
```

## 주의사항

- **conda 환경**: `snakemake`가 PATH에 있어야 하므로 실행 전 conda 환경 활성화 필수
- **병렬 실행 없음**: Snakemake 자체가 멀티코어를 사용하므로 프로젝트 간 병렬화는 리소스 경합 발생
- **실행 위치**: 파이프라인 루트(`/data_3tb/shared/rna-seq-pipeline`)에서 실행 권장
- **nohup 사용**: 장시간 실행 시 터미널 종료에 대비해 nohup 또는 tmux 사용 권장

```bash
# 백그라운드 실행 예시
nohup python scripts/batch_runner.py config/batch/my_batch.yaml --cores 12 \
  > batch_logs/nohup.out 2>&1 &
echo "PID: $!"
```
