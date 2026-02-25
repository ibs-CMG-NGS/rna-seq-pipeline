#!/bin/bash
# Ollama Setup Script for NGS Pipelines
# Current: RNA-seq Pipeline (testing & standardization)
# Future: Extensible to ATAC-seq, WGS, etc.
# 
# Usage: bash scripts/setup_ollama.sh

set -e  # Exit on error

echo "=================================================="
echo "🤖 Ollama Setup for NGS Pipelines"
echo "   Phase 1: RNA-seq Pipeline"
echo "=================================================="
echo ""

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Functions
print_success() {
    echo -e "${GREEN}✅ $1${NC}"
}

print_error() {
    echo -e "${RED}❌ $1${NC}"
}

print_info() {
    echo -e "${YELLOW}ℹ️  $1${NC}"
}

# Check system requirements
echo "Step 1: Checking system requirements..."
echo ""

# Check OS
if [[ ! -f /etc/os-release ]]; then
    print_error "Cannot detect OS version"
    exit 1
fi

source /etc/os-release
print_success "OS: $NAME $VERSION"

# Check RAM
TOTAL_RAM=$(free -h | awk '/^Mem:/ {print $2}')
RAM_GB=$(free -g | awk '/^Mem:/ {print $2}')
print_info "Available RAM: $TOTAL_RAM"

if [[ $RAM_GB -lt 8 ]]; then
    print_error "Insufficient RAM. Need at least 8GB, have ${RAM_GB}GB"
    exit 1
fi
print_success "RAM check passed"

# Check disk space
FREE_SPACE=$(df -h /home | awk 'NR==2 {print $4}')
FREE_SPACE_GB=$(df -BG /home | awk 'NR==2 {print $4}' | sed 's/G//')
print_info "Free disk space: $FREE_SPACE"

if [[ $FREE_SPACE_GB -lt 10 ]]; then
    print_error "Insufficient disk space. Need at least 10GB, have ${FREE_SPACE_GB}GB"
    exit 1
fi
print_success "Disk space check passed"

echo ""

# Check if Ollama is already installed
echo "Step 2: Checking Ollama installation..."
echo ""

if command -v ollama &> /dev/null; then
    OLLAMA_VERSION=$(ollama --version 2>&1 | head -1)
    print_info "Ollama already installed: $OLLAMA_VERSION"
    
    read -p "Reinstall Ollama? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        print_info "Skipping Ollama installation"
        SKIP_INSTALL=true
    fi
fi

if [[ ! $SKIP_INSTALL ]]; then
    print_info "Installing Ollama..."
    
    # Download and install
    curl -fsSL https://ollama.com/install.sh | sh
    
    if command -v ollama &> /dev/null; then
        print_success "Ollama installed successfully"
        ollama --version
    else
        print_error "Ollama installation failed"
        exit 1
    fi
fi

echo ""

# Start Ollama service
echo "Step 3: Starting Ollama service..."
echo ""

# Check if systemd service exists
if systemctl list-unit-files | grep -q ollama.service; then
    print_info "Starting Ollama systemd service..."
    sudo systemctl start ollama
    sudo systemctl enable ollama
    sleep 2
    
    if systemctl is-active --quiet ollama; then
        print_success "Ollama service is running"
        
        # Check if custom model path is set
        OLLAMA_MODELS_PATH=$(systemctl show ollama -p Environment | grep OLLAMA_MODELS || echo "")
        if [[ -z "$OLLAMA_MODELS_PATH" ]]; then
            print_info "Consider setting custom model path for better disk management:"
            echo "    sudo mkdir -p /data_3tb/shared/ollama-models"
            echo "    sudo systemctl edit ollama"
            echo "    # Add: Environment=\"OLLAMA_MODELS=/data_3tb/shared/ollama-models\""
        fi
    else
        print_error "Ollama service failed to start"
        sudo journalctl -u ollama -n 20 --no-pager
        exit 1
    fi
else
    print_info "Systemd service not found. Starting Ollama manually..."
    print_info "Run 'ollama serve' in a separate terminal"
    print_info "Or set up systemd service manually (see docs/developer/PHASE5_OLLAMA_SETUP.md)"
fi

# Verify Ollama is responding
print_info "Checking Ollama API..."
sleep 1

if curl -s http://localhost:11434/api/version > /dev/null; then
    print_success "Ollama API is responding"
else
    print_error "Cannot connect to Ollama API"
    print_info "Try starting manually: ollama serve"
    exit 1
fi

echo ""

# Pull Llama 3.1 model
echo "Step 4: Downloading Llama 3.1 model..."
echo ""

if ollama list | grep -q llama3.1:8b; then
    print_success "llama3.1:8b already downloaded"
else
    print_info "Downloading llama3.1:8b (~4.7GB)..."
    print_info "This may take several minutes..."
    
    ollama pull llama3.1:8b
    
    if ollama list | grep -q llama3.1:8b; then
        print_success "Model downloaded successfully"
    else
        print_error "Model download failed"
        exit 1
    fi
fi

echo ""

# Test Ollama
echo "Step 5: Testing Ollama..."
echo ""

print_info "Running simple test query..."
RESPONSE=$(ollama run llama3.1:8b "Say 'Hello from Ollama' and nothing else" --verbose=false 2>&1 | head -1)

if [[ -n "$RESPONSE" ]]; then
    print_success "Ollama test passed"
    print_info "Response: $RESPONSE"
else
    print_error "Ollama test failed"
    exit 1
fi

echo ""

# Summary
echo "=================================================="
echo "✅ Ollama Setup Complete!"
echo "=================================================="
echo ""
echo "Installation Summary:"
echo "  • Ollama version: $(ollama --version 2>&1 | head -1)"
echo "  • Model: llama3.1:8b"
echo "  • API endpoint: http://localhost:11434"
echo "  • Current pipeline: RNA-seq"
echo ""
echo "Next Steps:"
echo "  1. Test LLM agent (RNA-seq):"
echo "     cd /data_3tb/shared/rna-seq-pipeline"
echo "     conda activate rna-seq-pipeline"
echo "     python scripts/standardization/llm_agent.py \\"
echo "       --project-summary /data_3tb/shared/output/mouse-chd8/project_summary.json \\"
echo "       --rnaseq-output /data_3tb/shared/output/mouse-chd8"
echo ""
echo "  2. Try natural language queries:"
echo "     - 'mouse-chd8 프로젝트 QC 상태 보여줘'"
echo "     - '경로 검증해줘'"
echo "     - 'DE 분석 준비해줘'"
echo ""
echo "  3. View documentation:"
echo "     cat docs/developer/PHASE5_OLLAMA_SETUP.md"
echo ""
echo "Future Expansion:"
echo "  • ATAC-seq pipeline integration"
echo "  • WGS pipeline integration"
echo "  • Unified NGS agent (multi-pipeline)"
echo ""
echo "=================================================="
