# LLM Agent Security Guide

## 보안 위협 분석 및 대책

### 🔴 주요 보안 위협

#### 1. Command Injection
**위험**: LLM이 악의적인 명령어를 삽입할 수 있음

```python
# 공격 예시
project_id = "mouse-chd8; rm -rf /"
# 결과: subprocess.run에서 실행 → 파일 삭제 시도
```

**대책**: ✅ 구현됨
- Input validation with regex (`^[a-zA-Z0-9_-]+$`)
- Whitelist 기반 검증
- `scripts/utils/security.py` 사용

#### 2. Path Traversal
**위험**: 임의 경로 접근 시도

```python
# 공격 예시
project_id = "../../etc/passwd"
sample_id = "../../../root/.ssh/id_rsa"
```

**대책**: ✅ 구현됨
- Path traversal 패턴 차단 (`..`, `/`, `\\`)
- Allowed base directories 검증
- `validate_path()` 함수 사용

#### 3. Resource Exhaustion
**위험**: 무한 루프, 메모리 폭탄, Fork 폭탄

**대책**: ⚠️ 권장사항 (아직 미구현)
- systemd service로 실행 + resource limits
- Docker/Podman container 사용
- cgroups로 CPU/메모리 제한

#### 4. Arbitrary File Access
**위험**: 민감한 파일 읽기/쓰기

**대책**: ✅ 부분 구현
- Allowed directories whitelist
- Read-only 모드 권장
- Subprocess working directory 제한

## 구현된 보안 기능

### 1. Input Validation (`scripts/utils/security.py`)

```python
from scripts.utils.security import validate_project_id, SecurityError

try:
    project_id = validate_project_id(user_input)
except SecurityError as e:
    return {"status": "error", "message": str(e)}
```

**검증 항목**:
- ✅ Alphanumeric + hyphen + underscore only
- ✅ No path traversal (`..`, `/`, `\\`)
- ✅ Length limit (100 chars)
- ✅ No shell metacharacters (`;`, `|`, `&`, `$`, etc.)

### 2. Path Validation

```python
from scripts.utils.security import validate_path

allowed_dirs = [
    "/data_3tb/shared/output",
    "/data_3tb/shared/rna-seq-pipeline",
    "/home/ygkim/ngs-pipeline"
]

safe_path = validate_path(user_path, allowed_dirs)
```

### 3. Command Sanitization

```python
from scripts.utils.security import sanitize_command_args

args = ["--project-id", project_id, "--skip-de"]
safe_args = sanitize_command_args(args)
```

**차단되는 패턴**:
- `;` - Command chaining
- `|` - Pipes
- `&` - Background execution
- `$` - Variable expansion
- `` ` `` - Command substitution
- `()` - Subshells
- `>`, `<` - Redirection

## 권장 배포 구성

### Option 1: Systemd Service (추천)

```ini
# /etc/systemd/system/rna-seq-agent.service
[Unit]
Description=RNA-seq LLM Agent
After=network.target ollama.service

[Service]
Type=simple
User=ygkim
Group=lab
WorkingDirectory=/data_3tb/shared/rna-seq-pipeline

# Resource Limits
CPUQuota=400%              # Max 4 cores
MemoryLimit=32G            # Max 32GB RAM
TasksMax=100               # Max 100 processes
LimitNOFILE=1024           # Max open files

# Security
PrivateTmp=yes             # Isolated /tmp
ProtectSystem=strict       # Read-only /usr, /boot, /efi
ProtectHome=yes            # No access to /home (except WorkingDirectory)
ReadWritePaths=/data_3tb/shared  # Only write to data dir
NoNewPrivileges=yes        # Can't escalate privileges

# Execution
ExecStart=/opt/conda/envs/rna-seq-pipeline/bin/python \
    /data_3tb/shared/rna-seq-pipeline/scripts/standardization/llm_agent.py \
    --project-summary /data_3tb/shared/output/mouse-chd8/project_summary.json \
    --rnaseq-output /data_3tb/shared/output/mouse-chd8

Restart=on-failure
RestartSec=10s

[Install]
WantedBy=multi-user.target
```

**설치**:
```bash
sudo cp rna-seq-agent.service /etc/systemd/system/
sudo systemctl daemon-reload
sudo systemctl enable rna-seq-agent
sudo systemctl start rna-seq-agent
```

### Option 2: Docker Container (최고 보안)

```dockerfile
# Dockerfile
FROM condaforge/mambaforge:latest

# Install dependencies
COPY environment.yaml /tmp/
RUN mamba env create -f /tmp/environment.yaml && mamba clean -afy

# Copy code (read-only)
COPY --chown=nobody:nogroup . /app
WORKDIR /app

# Run as non-root
USER nobody

# Resource limits set via docker run
ENTRYPOINT ["conda", "run", "-n", "rna-seq-pipeline", "python", \
            "scripts/standardization/llm_agent.py"]
```

**실행**:
```bash
docker build -t rna-seq-agent .

docker run -it --rm \
    --cpus=4 \
    --memory=32g \
    --pids-limit=100 \
    --read-only \
    --tmpfs /tmp \
    -v /data_3tb/shared/output:/data:ro \
    -v /data_3tb/shared/rna-seq-pipeline:/app:ro \
    rna-seq-agent \
    --project-summary /data/mouse-chd8/project_summary.json \
    --rnaseq-output /data/mouse-chd8
```

### Option 3: Podman (rootless, 더 안전)

```bash
# Same as Docker but runs without root
podman run -it --rm \
    --cpus=4 \
    --memory=32g \
    --pids-limit=100 \
    --read-only \
    --tmpfs /tmp \
    -v /data_3tb/shared/output:/data:ro \
    -v /data_3tb/shared/rna-seq-pipeline:/app:ro \
    rna-seq-agent
```

## 추가 권장사항

### 1. Ollama Security

```yaml
# Ollama는 localhost:11434에서만 접근
OLLAMA_HOST=127.0.0.1:11434
OLLAMA_MODELS=/var/lib/ollama/models  # 모델은 로컬 저장

# 외부 네트워크 차단
sudo ufw deny from any to any port 11434
sudo ufw allow from 127.0.0.1 to 127.0.0.1 port 11434
```

### 2. Audit Logging

모든 LLM agent 실행을 로깅:

```python
import logging

logging.basicConfig(
    filename='/var/log/rna-seq-agent.log',
    level=logging.INFO,
    format='%(asctime)s - %(user)s - %(action)s - %(args)s'
)

# Tool 실행 전
logging.info(f"User: {user}, Tool: {tool_name}, Args: {arguments}")
```

### 3. Rate Limiting

DDoS/abuse 방지:

```python
from functools import wraps
import time

def rate_limit(max_calls=10, period=60):
    """Max 10 calls per minute"""
    calls = []
    
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            now = time.time()
            calls[:] = [c for c in calls if c > now - period]
            
            if len(calls) >= max_calls:
                raise Exception("Rate limit exceeded")
            
            calls.append(now)
            return func(*args, **kwargs)
        
        return wrapper
    return decorator

@rate_limit(max_calls=10, period=60)
def execute_tool(tool_name, args):
    ...
```

### 4. 사용자 권한 분리

```bash
# 전용 사용자 생성
sudo useradd -r -s /bin/false rna-seq-agent

# 최소 권한만 부여
sudo chown -R rna-seq-agent:lab /data_3tb/shared/rna-seq-pipeline
sudo chmod 550 /data_3tb/shared/rna-seq-pipeline  # Read + execute only
```

## 테스트

### Security Test Suite

```bash
# 실행
python scripts/utils/security.py

# 예상 출력:
# ✅ Valid project_id: mouse-chd8
# ✅ Caught: SecurityError("Invalid project_id")
# ✅ Caught: SecurityError("Path traversal detected")
# ✅ All security tests passed
```

### Penetration Testing

```python
# Test command injection
test_inputs = [
    "mouse-chd8; rm -rf /",
    "mouse-chd8 && cat /etc/passwd",
    "../../etc/shadow",
    "mouse-chd8`whoami`",
    "$(cat /etc/passwd)",
]

for inp in test_inputs:
    try:
        validate_project_id(inp)
        print(f"❌ FAILED: {inp} should be blocked!")
    except SecurityError:
        print(f"✅ PASSED: {inp} correctly blocked")
```

## 현재 상태 요약

### ✅ 구현됨
- Input validation (project_id, sample_id)
- Path traversal 방지
- Shell metacharacter 차단
- Allowed directories whitelist

### ⚠️ 권장사항 (미구현)
- systemd resource limits
- Docker/Podman containerization
- Audit logging
- Rate limiting
- User separation

### 🔴 TODO (Phase 4+)
- Automated security testing
- Penetration testing
- Security audit
- Documentation review

## 결론

**현재 보안 수준**: ⭐⭐⭐☆☆ (중간)

**프로덕션 배포 전 필수**:
1. systemd service 설정 (resource limits)
2. Audit logging 활성화
3. 전용 사용자 생성
4. Rate limiting 구현

**최고 보안을 위해**:
- Docker/Podman container 사용
- Read-only 파일시스템
- Network isolation (localhost only for Ollama)
- Regular security audits
