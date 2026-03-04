# LLM Agent 빠른 시작 가이드

RNA-seq 파이프라인을 자연어로 제어하는 AI agent 사용법입니다.  
**FASTQ 감지 → 프로젝트 설정 → 검증 → 실행 → 모니터링** 전 과정을 대화로 처리할 수 있습니다.

---

## 1. 사전 준비 (최초 1회)

```bash
# Ollama 설치 (서버에 이미 설치된 경우 생략)
curl -fsSL https://ollama.com/install.sh | sh

# Ollama 서비스 시작
ollama serve &

# 권장 모델 다운로드 (GPU 서버 기준, ~20GB)
ollama pull qwen2.5:32b

# Python 패키지 설치
conda activate rna-seq-pipeline
pip install ollama
```

> **서버 공용 이용 시**: Ollama는 이미 실행 중입니다. `ollama list`로 모델을 확인하세요.

---

## 2. Agent 실행

```bash
cd /home/ygkim/ngs-pipeline/rna-seq-pipeline

# 기본 (파이프라인 실행 목적)
python scripts/standardization/llm_agent.py \
  --interactive \
  --model qwen2.5:32b

# 기존 프로젝트 QC 분석 목적
python scripts/standardization/llm_agent.py \
  --project-summary /data_3tb/shared/output/mouse-chd8/project_summary.json \
  --rnaseq-output /data_3tb/shared/output/mouse-chd8 \
  --interactive \
  --model qwen2.5:32b
```

---

## 3. 전체 워크플로우 예시 (처음부터 끝까지)

### Step 1 — FASTQ 파일 확인
```
You: /data_3tb/shared/chd8-rna-seq-raw-data/fastq/ 에서 FASTQ 파일 찾아줘

Agent: 76개 파일, 38개 샘플 발견 (149 GB, paired-end)
       파일 패턴: Chd8_{tissue}_{age}_{sex}_{replicate}_1/2.fq.gz
```

### Step 2 — 프로젝트 설정 생성
```
You: 새 프로젝트 만들어줘. ID는 mouse-chd8-2026,
     결과는 /data_3tb/shared/output/mouse-chd8-2026/ 에 저장

Agent: config/projects/mouse-chd8-2026.yaml 생성 완료 ✅
```

### Step 3 — 리소스 확인
```
You: 리소스 얼마나 필요해?

Agent: 입력 149GB → 예상 출력 746GB
       소요 시간: 5-10시간 (8 cores)
       RAM: ~30GB (STAR), 디스크 여유 2145GB ✅
```

### Step 4 — 샘플 시트 생성
```
You: 샘플 시트 만들어줘

Agent: 38개 샘플 자동 분류 완료
       wildtype: 20개, heterozygous: 18개
       → config/samples/mouse-chd8-2026.tsv 저장
```

### Step 5 — 입력 데이터 검증
```
You: config/projects/mouse-chd8-2026.yaml 로 데이터 검증해줘

Agent: ✅ FASTQ 파일 존재 확인
       ✅ 디스크 여유 공간 충분
       ⚠️  genome 경로 확인 필요: /data_3tb/shared/genome/star_index/
```

### Step 6 — Dry-run
```
You: dry-run 해줘

Agent: 270개 jobs, 11개 rule 확인
       cutadapt: 38 / fastqc_raw: 76 / star_align: 38 / ...
       문제 없음 ✅ 실제 실행 준비 완료
```

### Step 7 — 파이프라인 실행
```
You: 실제로 실행해줘, 16 cores로

Agent: 파이프라인 시작 🚀
       snakemake --configfile ... --cores 16
       진행 상황은 "파이프라인 상태 보여줘" 로 확인하세요
```

### Step 8 — 모니터링
```
You: 파이프라인 상태 보여줘

Agent: 진행률: 45% (cutadapt 완료, star_align 진행 중)
       완료: cutadapt(38/38), fastqc(76/76)
       진행 중: star_align (17/38)
```

### Step 9 — 완료 후 QC 확인
```
You: QC 결과 어때?

Agent: 38개 샘플 전부 통과 (100%)
       평균 맵핑률: 92.4%, 평균 assignment: 69.9%
       DE 분석 준비 완료 ✅
```

---

## 4. QC 분석 쿼리 예시

```
# 전체 현황
You: QC 상태 보여줘
You: 몇 개 샘플 통과했어?

# 조건 비교
You: wildtype vs heterozygous 비교해줘
You: HPC에서만 비교해줘
You: Male과 Female 비교

# 다차원 분석
You: 어떤 그룹(축)이 있어?
You: HPC wildtype male 샘플 목록 보여줘

# 실패 샘플
You: QC 실패한 샘플 있어?
You: Chd8_HPC_10M_W_1 샘플 상세 정보

# DE 분석 연동
You: DE 분석 시작해줘
```

---

## 5. 명령줄 옵션

| 옵션 | 설명 | 예시 |
|------|------|------|
| `--interactive` | 대화형 모드 (권장) | — |
| `--message "..."` | 단일 쿼리 모드 | `--message "QC 상태 보여줘"` |
| `--model` | 모델 선택 | `--model qwen2.5:32b` |
| `--project-summary` | 기존 프로젝트 요약 경로 | `--project-summary .../project_summary.json` |
| `--rnaseq-output` | 파이프라인 결과 디렉토리 | `--rnaseq-output /data_3tb/.../output` |
| `--ollama-host` | 원격 Ollama 서버 주소 | `--ollama-host http://gpu-server:11434` |

---

## 6. 모델 선택 가이드

| 모델 | VRAM/RAM | 속도 | 품질 | Tool Calling | 권장 용도 |
|------|----------|------|------|--------------|-----------|
| **qwen2.5:32b** | 20GB VRAM | ⚡⚡ | ⭐⭐⭐⭐⭐ | Native (GPU) | **기본 권장** |
| llama3.1:8b | 8GB RAM | ⚡⚡⚡ | ⭐⭐⭐⭐ | Native | GPU 없을 때 |
| mistral:7b | 8GB RAM | ⚡⚡⚡⚡ | ⭐⭐⭐ | Text pattern | 빠른 상태 확인 |
| llama3.1:70b | 42GB RAM | ⚡ | ⭐⭐⭐⭐⭐ | Native | 고성능 서버 |

---

## 7. 사용 가능한 도구 (17개)

**QC 분석 (Phase 6)**
- `get_project_status` — 전체 QC 현황
- `compare_conditions` — 조건 간 비교
- `get_failed_samples` — 실패 샘플 목록
- `get_sample_details` — 샘플 상세 정보
- `list_conditions` — 실험 조건 목록
- `start_de_analysis` — DE/GO 분석 시작

**다차원 분석 (Phase 6.1)**
- `get_sample_axes` — 실험 축 목록 (tissue, sex, genotype, ...)
- `compare_by_axis` — 임의 축으로 비교 + 필터링
- `filter_samples` — 다중 조건 필터링

**파이프라인 실행 (Phase 8A)**
- `create_project_config` — 프로젝트 config 생성
- `detect_fastq_files` — FASTQ 파일 감지
- `validate_input_data` — 사전 검증
- `run_pipeline` — Snakemake 실행 (dry-run / 실제)

**모니터링/리소스 (Phase 8B)**
- `monitor_pipeline` — 실행 상태 확인
- `create_sample_sheet` — 샘플 시트 자동 생성
- `estimate_resources` — 리소스 예측 (시간/디스크/RAM)

---

## 8. 문제 해결

**Ollama 연결 안 됨**
```bash
ollama serve          # 서비스 시작
ollama list           # 모델 목록 확인
```

**모델 없음**
```bash
ollama pull qwen2.5:32b    # 권장 모델 (~20GB)
ollama pull llama3.1:8b    # 경량 대안 (~5GB)
```

**응답이 느림**
```bash
nvidia-smi            # GPU 인식 확인
# GPU 없으면 경량 모델 사용: --model llama3.1:8b
```

**config 파일 경로 오류**  
프로젝트 루트 기준 상대경로 또는 절대경로 모두 사용 가능합니다.

---

## 상세 문서

- Ollama 설치/설정: `docs/user/OLLAMA_SETUP.md`
- 파이프라인 전체 가이드: `docs/user/PIPELINE_GUIDE.md`
- 개발자 기술 문서: `docs/developer/PHASE6_LLM_INTEGRATION.md`
- Phase 8 파이프라인 실행: `docs/developer/PHASE8_PIPELINE_EXECUTION.md`
