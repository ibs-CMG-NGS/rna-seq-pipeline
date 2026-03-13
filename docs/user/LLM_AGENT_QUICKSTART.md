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

### 방법 A — VS Code Task (권장)
VS Code에서 `Ctrl+Shift+B` → **RNA-seq Agent (interactive)** 선택

기존 프로젝트 이어서 시작: **RNA-seq Agent (with existing project)** → config 파일 경로 입력
세션 초기화 후 시작: **RNA-seq Agent (reset session)**

### 방법 B — 터미널
```bash
cd /data_3tb/shared/rna-seq-pipeline
conda activate rna-seq-pipeline

# 기본 (새 프로젝트 또는 이어서)
python scripts/standardization/llm_agent.py \
  --interactive \
  --model qwen3:14b

# 기존 프로젝트 config 지정
python scripts/standardization/llm_agent.py \
  --interactive \
  --model qwen3:14b \
  --config config/projects/mouse-chd8.yaml

# 세션 초기화 후 시작
python scripts/standardization/llm_agent.py \
  --interactive \
  --model qwen3:14b \
  --reset-session
```

---

## 3. 전체 워크플로우 예시 (처음부터 끝까지)

### 권장 2단계 워크플로우 (2026-03-11 이후)

**Stage 1 — 한번에 설정** (FASTQ 감지 + 샘플시트 생성 + 라이브러리 검증 통합)
```
You: 프로젝트 설정해줘

Agent: [FASTQ] 76개 파일, 38개 샘플 감지 완료
       [샘플시트] 38개 샘플 등록 (wildtype: 20개, heterozygous: 18개)
       [라이브러리 검증] 샘플: Chd8_HPC_10M_W_1, strandedness: reverse (2),
                        엑손 비율: 71.3%, 판정: mRNA-seq ✅
       → strandedness config 자동 업데이트: 0 → 2
```

**Stage 2 — 파이프라인 실행**
```
You: dry-run 해줘

Agent: 270개 jobs, 11개 rule 확인. 문제 없음 ✅

You: 실제로 실행해줘, 16 cores로

Agent: 파이프라인 시작. 진행 상황은 "파이프라인 상태 보여줘" 로 확인
```

---

### 단계별 워크플로우 (상세 제어)

### Step 1 — FASTQ 파일 확인
```
You: /data_3tb/shared/chd8-rna-seq-raw-data/fastq/ 에서 FASTQ 파일 찾아줘

Agent: 76개 파일, 38개 샘플 발견 (149 GB, paired-end)
```

### Step 2 — 프로젝트 설정 생성
```
You: 새 프로젝트 만들어줘. ID는 mouse-chd8-2026,
     결과는 /data_3tb/shared/output/ 에 저장

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

### Step 5 — 라이브러리 검증 (파이프라인 실행 전 필수)
```
You: 라이브러리 검증해줘

Agent: strandedness: reverse (2), 엑손 비율: 71.3%, 판정: mRNA-seq ✅
       config 자동 업데이트 완료
```

### Step 6 — 입력 데이터 검증
```
You: 데이터 검증해줘

Agent: ✅ FASTQ 파일 존재 확인
       ✅ 디스크 여유 공간 충분
       ⚠️  genome 경로 확인 필요: /data_3tb/shared/genome/star_index/
```

### Step 7 — Dry-run
```
You: dry-run 해줘

Agent: 270개 jobs, 11개 rule 확인. 문제 없음 ✅
```

### Step 8 — 파이프라인 실행
```
You: 실제로 실행해줘, 16 cores로

Agent: 파이프라인 시작. 진행 상황은 "파이프라인 상태 보여줘" 로 확인하세요
```

### Step 9 — 모니터링
```
You: 파이프라인 상태 보여줘

Agent: 진행률: 45% (cutadapt 완료, star_align 진행 중)
       완료: cutadapt(38/38), fastqc(76/76)
       진행 중: star_align (17/38)
```

### Step 10 — 완료 후 QC 확인
```
You: QC 결과 어때?

Agent: 38개 샘플 전부 통과 (100%)
       평균 맵핑률: 92.4%, 평균 assignment: 69.9%
       DE 분석 준비 완료 ✅
```

### Step 11 — DE-GO 파이프라인 준비 (bridge)
```
You: DE 분석 연결해줘

Agent: counts matrix / metadata / DE config 복사 완료
       DE config 조건명 자동 검증... 불일치 없음 ✅
```

### Step 12 — DE-GO 파이프라인 실행
```
You: DE 분석 실행해줘

Agent: [dry-run] 172 jobs, 16개 rule 확인
       run_pairwise_de(4), go_enrichment(36), kegg_enrichment(12) 외 다수
       실제로 실행할까요?

You: 응

Agent: DE 파이프라인 백그라운드 시작 (PID: 12345)
       진행 상황은 "DE 진행상황 알려줘" 로 확인하세요
```

### Step 13 — DE-GO 모니터링 및 결과 확인
```
You: DE 진행상황 알려줘

Agent: 진행률: 62% (38/61 jobs) — 현재: go_enrichment
       완료: run_pairwise_de(4), generate_pairwise_qc(4) 외 다수

You: DE 결과 요약해줘

Agent: [D1_vs_Control] DEG 823개 (up: 412, down: 411)
                       상위 GO: oxidative stress response, apoptosis ...
       [24h_vs_Control] DEG 1204개 ...
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
| `--model` | 모델 선택 | `--model qwen3:14b` |
| `--config` | 기존 프로젝트 config 경로 (세션 복원) | `--config config/projects/mouse-chd8.yaml` |
| `--reset-session` | 이전 세션 초기화 후 시작 | — |
| `--project-summary` | 기존 프로젝트 요약 경로 | `--project-summary .../project_summary.json` |
| `--rnaseq-output` | 파이프라인 결과 디렉토리 | `--rnaseq-output /data_3tb/.../output` |
| `--ollama-host` | 원격 Ollama 서버 주소 | `--ollama-host http://gpu-server:11434` |

> **세션 지속성**: Agent는 대화 종료 후에도 현재 프로젝트 정보(config 경로, data_dir 등)를 `.agent_session.json`에 저장합니다. 다음 실행 시 자동으로 복원됩니다.

---

## 6. 모델 선택 가이드

| 모델 | VRAM/RAM | 속도 | 품질 | Tool Calling | 권장 용도 |
|------|----------|------|------|--------------|-----------|
| **qwen3:14b** | ~10GB RAM | ⚡⚡⚡ | ⭐⭐⭐⭐⭐ | Native | **현재 기본 권장** |
| qwen2.5:32b | 20GB VRAM | ⚡⚡ | ⭐⭐⭐⭐⭐ | Native (GPU) | GPU 서버에서 고품질 |
| llama3.1:8b | 8GB RAM | ⚡⚡⚡ | ⭐⭐⭐⭐ | Native | 경량 대안 |
| mistral:7b | 8GB RAM | ⚡⚡⚡⚡ | ⭐⭐⭐ | Text pattern | 빠른 상태 확인 |

---

## 7. 사용 가능한 도구 (28개)

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

**라이브러리 검증 (Phase 9)**
- `validate_library_type` — RSeQC로 strandedness 자동 감지 + mRNA-seq 여부 확인, config 자동 업데이트
- `setup_and_validate` — FASTQ 감지 + 샘플시트 + 라이브러리 검증 한번에 실행 (새 프로젝트 권장)

**DE-GO 파이프라인 연결 (Phase 10)**
- `run_bridge` — RNA-seq 결과를 DE 파이프라인으로 연결 (counts/metadata/config 복사)
- `validate_de_config_conditions` — DE config 조건명 검증 (metadata와 불일치 자동 감지)
- `apply_de_config_corrections` — 조건명 불일치 자동 수정
- `check_analysis_readiness` — DE 실행 전 5-point 사전 점검

**DE-GO 파이프라인 실행 (Phase 10)**
- `run_de_pipeline` — DE-GO Snakemake 파이프라인 실행 (dry-run / 실제 / 백그라운드)
- `monitor_de_pipeline` — 실행 중인 DE 파이프라인 진행상황 확인

**DE 결과 해석 (Phase 10)**
- `read_de_results_summary` — DEG 수, 상위 유전자, GO 결과 요약 (인지 도구)
- `read_counts` — 유전자별 발현량 조회

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
