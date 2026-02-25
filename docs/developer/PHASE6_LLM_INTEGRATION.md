# Phase 6: LLM Agent Integration Guide

## 🎯 목표

명령어 기반 쿼리를 **자연어 대화**로 전환하여 직관적인 파이프라인 관리 가능

## 📋 사전 준비

### 1. LLM API Key 준비

**Option A: OpenAI GPT-4** (권장)
```bash
# API key 발급: https://platform.openai.com/api-keys
export OPENAI_API_KEY="sk-..."

# Python 패키지 설치
pip install openai
```

**Option B: Anthropic Claude**
```bash
# API key 발급: https://console.anthropic.com/
export ANTHROPIC_API_KEY="sk-ant-..."

# Python 패키지 설치
pip install anthropic
```

### 2. 의존성 설치

```bash
conda activate rna-seq-pipeline
pip install openai  # 또는 anthropic
```

## 💬 사용법

### Interactive Mode (대화형)

```bash
python scripts/standardization/llm_agent.py \
    --project-summary /data/mouse-chd8/project_summary.json \
    --rnaseq-output /data/mouse-chd8 \
    --de-pipeline ../RNA-Seq_DE_GO_analysis \
    --interactive
```

**대화 예시**:

```
🤖 RNA-seq Pipeline Agent
Project: mouse-chd8
Type 'exit' or 'quit' to end session

You: What's the QC status?

Agent: Your project has 38 samples with a 100% QC pass rate! 
All samples passed quality control successfully. 

The samples are divided into:
- wildtype: 20 samples
- heterozygous: 18 samples

Would you like me to compare the QC metrics between these conditions?

You: Yes, compare them

Agent: Here's the QC comparison between wildtype and heterozygous samples:

**Mapping Rate:**
- Wildtype: 92.39% (average)
- Heterozygous: 92.49% (average)
- Difference: Negligible (0.1%)

**Assignment Rate:**
- Wildtype: 73.37% (average)
- Heterozygous: 73.51% (average)
- Difference: Negligible (0.14%)

Both groups show excellent and nearly identical quality metrics. 
There's no indication of systematic technical bias between conditions.
You can confidently proceed with differential expression analysis.

You: Start DE analysis

Agent: I'll start the differential expression and GO enrichment analysis.

First, let me prepare:
1. ✅ Copying counts matrix
2. ✅ Generating metadata from sample sheet
3. ✅ Creating config file (config_mouse-chd8.yml)
   - Species: mouse (auto-detected)
   - Comparison: heterozygous_vs_wildtype

Ready to start? This will:
- Run DESeq2 differential expression analysis
- Perform GO and KEGG enrichment analysis
- Generate publication-ready plots

Shall I proceed? (This will take ~30-60 minutes)

You: Yes

Agent: Starting analysis now... 
[Bridge script executing]
✅ Analysis started successfully!

You can monitor progress at:
  ../RNA-Seq_DE_GO_analysis/output/mouse-chd8/

I'll notify you when it's complete.
```

### Single Query Mode

```bash
# 단일 질문
python scripts/standardization/llm_agent.py \
    --project-summary /data/mouse-chd8/project_summary.json \
    --rnaseq-output /data/mouse-chd8 \
    --message "Which samples failed QC?"

# 출력: All 38 samples passed QC. No failures detected.
```

## 🔧 고급 설정

### 커스텀 프롬프트

`llm_agent.py`의 `_build_system_prompt()` 메서드를 수정하여 에이전트 행동 조정:

```python
def _build_system_prompt(self) -> str:
    return f"""You are a bioinformatics expert specializing in RNA-seq analysis.

Project: {self.project_id}
Your expertise: QC interpretation, statistical analysis, GO enrichment

Communication style:
- Be precise with numbers but explain their biological significance
- Suggest next steps proactively
- Flag potential issues early
- Use emojis sparingly for clarity (✅ ❌ ⚠️)

When user asks about QC:
1. Give clear pass/fail status
2. Explain any issues in biological context
3. Recommend remediation if needed
"""
```

### 추가 Tools 정의

새로운 기능 추가:

```python
{
    "name": "rerun_failed_samples",
    "description": "Rerun RNA-seq pipeline for failed samples",
    "parameters": {
        "type": "object",
        "properties": {
            "sample_ids": {
                "type": "array",
                "items": {"type": "string"},
                "description": "List of sample IDs to rerun"
            }
        },
        "required": ["sample_ids"]
    }
}
```

## 🎓 LLM Function Calling 작동 원리

### 1. User Query → LLM

```
User: "Compare wildtype vs heterozygous"
  ↓
LLM: "I need to use the 'compare_conditions' tool"
```

### 2. LLM → Tool Execution

```python
# LLM이 생성한 함수 호출
{
    "function_call": {
        "name": "compare_conditions",
        "arguments": "{}"
    }
}

# Python이 실제 실행
result = agent._execute_tool("compare_conditions", {})
# → subprocess로 agent_query.py 실행
```

### 3. Tool Result → LLM → Natural Language

```python
# Tool 결과 (JSON)
{
    "wildtype": {"mean_uniquely_mapped": 92.39, ...},
    "heterozygous": {"mean_uniquely_mapped": 92.49, ...}
}

# LLM이 자연어로 변환
"Both groups show excellent mapping rates (92.39% vs 92.49%). 
The small difference is not statistically or biologically significant..."
```

## 🚀 실전 워크플로우

### 시나리오: 새로운 RNA-seq 프로젝트 분석

```bash
# 1. RNA-seq 파이프라인 실행 (기존)
snakemake --cores 12

# 2. Agent 시작
python scripts/standardization/llm_agent.py \
    --project-summary results/project_summary.json \
    --rnaseq-output results/ \
    --de-pipeline ../RNA-Seq_DE_GO_analysis \
    --interactive
```

**대화 흐름**:

```
You: Summarize the project

Agent: [자동으로 get_project_status 호출]
       "This is a mouse CHD8 study with 38 samples..."

You: Any quality concerns?

Agent: [자동으로 get_failed_samples 호출]
       "All samples passed. Mean mapping 92.4%..."

You: Good. Start DE analysis with default settings

Agent: [자동으로 start_de_analysis 호출]
       [실제로 bridge_to_de_pipeline.py 실행]
       "Analysis started..."

You: When will it finish?

Agent: "Based on 38 samples, estimated 45-60 minutes.
       I'll check the output directory for completion..."
```

## ⚠️ 주의사항

### 1. API 비용

- GPT-4: ~$0.03 per 1K tokens
- 대화 세션당 예상: $0.10-0.50
- Function calling은 추가 토큰 소비

### 2. 보안

```bash
# API key를 환경 변수로 관리
export OPENAI_API_KEY="sk-..."

# 코드에 직접 넣지 말 것!
# ❌ api_key = "sk-proj-..."
```

### 3. Rate Limiting

- OpenAI: 분당 요청 수 제한
- 긴 대화는 rate limit에 주의

### 4. 에러 핸들링

```python
try:
    response = agent.chat(user_input)
except openai.error.RateLimitError:
    print("API rate limit reached. Wait 60 seconds...")
except openai.error.APIError as e:
    print(f"API error: {e}")
```

## 📊 Phase 6 체크리스트

- [ ] LLM API key 발급 및 설정
- [ ] `llm_agent.py` 실행 환경 구축
- [ ] 기본 대화 테스트 (status, compare)
- [ ] DE 분석 자동 실행 테스트
- [ ] 커스텀 프롬프트 조정
- [ ] 추가 tools 정의 (필요 시)
- [ ] 에러 핸들링 강화
- [ ] 사용자 문서 작성

## 🔮 미래 확장

### GitHub Copilot Chat Integration

VS Code에서 직접 대화:

```
User: @workspace RNA-seq QC 상태 보여줘

Copilot: [자동으로 llm_agent.py 실행]
         "38개 샘플 모두 통과..."
```

### Slack Bot Integration

```python
from slack_bolt import App

@app.message("qc status")
def handle_qc_status(message, say):
    agent = PipelineAgent(...)
    response = agent.chat(message['text'])
    say(response)
```

### Web UI

```python
# FastAPI backend
from fastapi import FastAPI

app = FastAPI()

@app.post("/chat")
async def chat(message: str):
    response = agent.chat(message)
    return {"response": response}
```

---

## 🎯 요약

**Phase 6는 선택적이지만 강력합니다**:

| 기능 | Phase 5까지 | Phase 6 추가 시 |
|:-----|:------------|:---------------|
| QC 확인 | `python agent_query.py --query status` | "QC 상태 보여줘" |
| 조건 비교 | `python agent_query.py --query compare` | "wildtype과 het 비교해줘" |
| DE 분석 | `python bridge_to_de_pipeline.py ...` | "DE 분석 시작" |
| 학습 곡선 | 명령어 문법 학습 필요 | 자연어만으로 가능 |
| 사용성 | 전문가용 | 누구나 사용 가능 |

**Phase 6 없이도 Phase 1-5는 완전히 작동합니다!**
