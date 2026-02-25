# LLM Agent - New DE/GO Bridge Tools

## Phase 3 완료: Bridge 관련 Tool 추가

LLM agent에 다음 4개의 새로운 tool이 추가되었습니다:

### 1. prepare_de_analysis
DE/GO 분석 준비 (counts 복사, metadata 생성, config 생성)

**사용 예시**:
```
User: "mouse-chd8 프로젝트를 DE 분석 준비해줘"
Agent: [calls prepare_de_analysis()]
Response: "✅ DE 분석 준비 완료. Counts, metadata, config 파일 생성됨"
```

### 2. check_bridge_config  
Bridge 설정 파일 존재 확인 및 생성

**사용 예시**:
```
User: "mouse-chd8의 bridge 설정 확인해줘"
Agent: [calls check_bridge_config(project_id="mouse-chd8")]
Response: "Config exists at /path/to/paths_mouse_chd8.yaml"
```

### 3. validate_paths
DE/GO 분석에 필요한 모든 경로 검증

**사용 예시**:
```
User: "DE 분석 실행 전에 경로 검증해줘"
Agent: [calls validate_paths(project_id="mouse-chd8")]
Response: "✅ rnaseq_output, ✅ counts_file, ✅ de_pipeline, ..."
```

### 4. start_de_analysis (업데이트)
기존 tool, 이제 자동으로 config 생성 후 실행

## 자연어 사용 예시

### 시나리오 1: 처음부터 DE 분석까지
```
User: "mouse-chd8 프로젝트로 DE/GO 분석하고 싶어"

Agent: 
1. [calls check_bridge_config()]
2. [calls validate_paths()]
3. [calls prepare_de_analysis()]
4. "준비 완료. DE 분석을 시작할까요?"

User: "응"

Agent:
5. [calls start_de_analysis(confirm=True)]
6. "DE/GO 분석이 시작되었습니다."
```

### 시나리오 2: 상태 확인 후 진행
```
User: "현재 QC 통과한 샘플로 DE 분석 가능해?"

Agent:
1. [calls get_project_status()]
2. "38개 샘플 중 38개 통과 (100%). DE 분석 가능합니다."

User: "그럼 바로 준비해줘"

Agent:
3. [calls prepare_de_analysis()]
4. "✅ 준비 완료. 비교 조건: heterozygous vs wildtype"
```

### 시나리오 3: 경로 문제 해결
```
User: "DE 분석 준비해줘"

Agent:
1. [calls validate_paths()]
2. "⚠️ 일부 경로 누락: sample_sheet_dir"

User: "자동으로 고쳐줘"

Agent:
3. [calls check_bridge_config(force_regenerate=True)]
4. [calls validate_paths()]
5. "✅ 모든 경로 검증 완료"
```

## Tool 함수 상세

### prepare_de_analysis()
```python
{
    "name": "prepare_de_analysis",
    "description": "Prepare DE/GO analysis by copying counts, generating metadata, and creating config",
    "parameters": {
        "project_id": "optional, uses current project if not specified"
    }
}
```

**Returns**:
```json
{
    "status": "success",
    "output": "... bridge script output ...",
    "message": "DE analysis preparation complete"
}
```

### check_bridge_config()
```python
{
    "name": "check_bridge_config",
    "parameters": {
        "project_id": "required",
        "force_regenerate": "optional, default False"
    }
}
```

**Returns**:
```json
{
    "status": "exists" | "created" | "error",
    "config_path": "/path/to/config.yaml",
    "message": "description"
}
```

### validate_paths()
```python
{
    "name": "validate_paths",
    "parameters": {
        "project_id": "required"
    }
}
```

**Returns**:
```json
{
    "status": "validated",
    "validation": {
        "rnaseq_output": true,
        "counts_file": true,
        "de_pipeline": true,
        ...
    },
    "all_valid": true,
    "missing_paths": []
}
```

## 테스트 방법

```bash
# 1. LLM agent 실행
cd /data_3tb/shared/rna-seq-pipeline
python scripts/standardization/llm_agent.py \
    --project-summary /data_3tb/shared/output/mouse-chd8/project_summary.json \
    --rnaseq-output /data_3tb/shared/output/mouse-chd8

# 2. 자연어로 테스트
> mouse-chd8 프로젝트 DE 분석 준비해줘
> 경로 검증해줘
> bridge 설정 확인해줘
```

## 기술 구현

### Tool Registration
`scripts/standardization/llm_agent.py`:
- `_define_tools()`: Tool 정의 추가
- `_execute_tool()`: Tool 실행 로직 추가

### 의존성
- `scripts/utils/auto_config.py`: Auto-config 유틸리티
- `scripts/bridge_to_de_pipeline.py`: Bridge 스크립트

### 통합 흐름
```
User Query → LLM → Function Call → Tool Execution → Result → LLM → Natural Response
```

## 다음 단계 (Phase 4)

E2E 테스트:
1. 자연어 입력
2. Agent가 자동으로 tool 선택
3. Path discovery
4. Config 생성
5. Bridge 실행
6. 결과 확인
