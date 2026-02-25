# Future: Unified NGS Agent Architecture

**Status**: Design Document  
**Created**: 2026-02-26  
**Current Phase**: RNA-seq agent testing (Phase 5-7)  
**Target**: Multi-pipeline integration (Phase 8+)

---

## Vision

**Goal**: Single LLM agent that can analyze any NGS data type with natural language interface

**User Experience:**
```
User: "여기 FASTQ 파일들 있어. 분석해줘."

Agent:
1. 🔍 Analyzing data characteristics...
   → Detected: RNA-seq, paired-end, 75bp reads
2. ✅ Selected pipeline: RNA-seq
3. 🧪 Running QC...
   → All samples passed quality filters
4. 📊 QC Summary: 38 samples, 100% pass rate
5. ❓ Proceed with downstream analysis? (DE/GO)
   → User: "응"
6. 🚀 Starting differential expression analysis...
```

---

## Architecture Overview

```
┌─────────────────────────────────────────────────────────────┐
│                    Unified NGS Agent                        │
│                  (Natural Language Interface)                │
└──────────────────────┬──────────────────────────────────────┘
                       │
                       ├─── Data Type Detector
                       │    ├─ Read length analysis
                       │    ├─ Fragment size distribution
                       │    ├─ Filename patterns
                       │    └─ User hints
                       │
                       ├─── Pipeline Router
                       │    ├─ RNA-seq → /data_3tb/shared/rna-seq-pipeline
                       │    ├─ ATAC-seq → /data_Xtb/atac-seq-pipeline
                       │    └─ WGS → /data_Ytb/wgs-pipeline
                       │
                       ├─── QC Manager
                       │    ├─ FastQC (all types)
                       │    ├─ Type-specific QC
                       │    └─ Unified reporting
                       │
                       └─── Downstream Coordinator
                            ├─ RNA-seq: DE analysis, GO enrichment
                            ├─ ATAC-seq: Peak calling, motif analysis
                            └─ WGS: Variant calling, annotation
```

---

## Current State (Phase 5-7)

### RNA-seq Agent (Active Development)

**Location**: `/data_3tb/shared/rna-seq-pipeline/scripts/standardization/llm_agent.py`

**Capabilities**:
- ✅ Get sample details and QC status
- ✅ Validate paths
- ✅ Prepare DE analysis
- ✅ Check bridge configuration
- ✅ Generate comparison groups
- ✅ Security validation

**Limitations**:
- Only handles RNA-seq data
- Assumes data is already processed (aligned, counted)
- No automatic pipeline execution from raw FASTQ

---

## Future Phases

### Phase 8: Data Type Detection Module

**Goal**: Automatically identify NGS data type from FASTQ files

**Components**:

```python
# detectors/data_type_detector.py

class DataTypeDetector:
    def detect(self, fastq_files):
        """
        Detect NGS data type from FASTQ characteristics.
        
        Returns:
            {
                'type': 'rna-seq',
                'confidence': 0.95,
                'evidence': {
                    'read_length': 75,
                    'paired_end': True,
                    'stranded': 'yes',
                    'filename_hint': True
                }
            }
        """
        features = self._extract_features(fastq_files)
        return self._classify(features)
```

**Detection Rules**:
```yaml
# configs/detection_rules.yaml

rna-seq:
  read_length: [50, 150]
  strandedness_check: true
  typical_patterns:
    - "*RNA*"
    - "*transcriptome*"
  fragment_size: [200, 500]  # for paired-end

atac-seq:
  read_length: [30, 100]
  fragment_size: bimodal  # nucleosome-free + mono-nucleosome
  typical_patterns:
    - "*ATAC*"
    - "*chromatin*"

wgs:
  read_length: [100, 250]
  coverage_expectation: genome-wide
  typical_patterns:
    - "*WGS*"
    - "*genome*"
    - "*DNA*"
```

---

### Phase 9: Pipeline Integration Layer

**Goal**: Unified interface for all pipelines

**Components**:

```python
# pipelines/base.py

class BasePipeline(ABC):
    """Abstract base class for all NGS pipelines."""
    
    @abstractmethod
    def validate_input(self, files):
        """Validate input files for this pipeline."""
        pass
    
    @abstractmethod
    def run_qc(self, files, output_dir):
        """Run QC analysis."""
        pass
    
    @abstractmethod
    def run_primary_analysis(self, files, config):
        """Run primary analysis (alignment, counting, etc.)."""
        pass
    
    @abstractmethod
    def run_downstream(self, primary_output, config):
        """Run downstream analysis."""
        pass
```

```python
# pipelines/rna_seq.py

class RNASeqPipeline(BasePipeline):
    def __init__(self):
        self.pipeline_dir = "/data_3tb/shared/rna-seq-pipeline"
        self.conda_env = "rna-seq-pipeline"
        self.workflow = "Snakefile"
    
    def run_qc(self, files, output_dir):
        # Run FastQC, MultiQC
        cmd = f"snakemake --cores 8 --until multiqc"
        return subprocess.run(cmd, shell=True)
    
    def run_downstream(self, primary_output, config):
        # Prepare DE/GO analysis
        bridge = BridgeToDEPipeline(primary_output)
        return bridge.execute()
```

---

### Phase 10: Unified Agent Implementation

**Location**: `/data_3tb/shared/ngs-agent/`

**Structure**:
```
ngs-agent/
├── main_agent.py              # Entry point
├── orchestrator.py            # Workflow orchestration
├── detectors/
│   ├── data_type_detector.py
│   └── quality_checker.py
├── pipelines/
│   ├── base.py
│   ├── rna_seq.py
│   ├── atac_seq.py
│   └── wgs.py
├── tools/
│   ├── file_validator.py
│   ├── qc_aggregator.py
│   └── report_generator.py
├── configs/
│   ├── pipeline_paths.yaml
│   ├── detection_rules.yaml
│   └── workflows.yaml
└── tests/
    ├── test_detector.py
    └── test_integration.py
```

**Main Agent**:
```python
# main_agent.py

class UnifiedNGSAgent:
    def __init__(self, ollama_endpoint="http://localhost:11434"):
        self.llm = OllamaClient(ollama_endpoint)
        self.detector = DataTypeDetector()
        self.pipelines = {
            'rna-seq': RNASeqPipeline(),
            'atac-seq': ATACSeqPipeline(),
            'wgs': WGSPipeline()
        }
        self.orchestrator = WorkflowOrchestrator()
    
    def process_request(self, user_query, data_files):
        """
        Process user's natural language request.
        
        Example:
            user_query = "이 데이터 분석해줘. QC 통과하면 downstream까지"
            data_files = ["sample1_R1.fastq.gz", "sample1_R2.fastq.gz"]
        """
        # 1. Parse user intent with LLM
        intent = self.llm.parse_intent(user_query)
        
        # 2. Detect data type
        data_type = self.detector.detect(data_files)
        self.llm.report(f"✅ Detected: {data_type['type']} "
                       f"(confidence: {data_type['confidence']:.0%})")
        
        # 3. Select pipeline
        pipeline = self.pipelines[data_type['type']]
        
        # 4. Run QC
        qc_result = pipeline.run_qc(data_files)
        qc_summary = self.llm.summarize_qc(qc_result)
        self.llm.report(f"📊 QC: {qc_summary}")
        
        # 5. Check if QC passed
        if not qc_result.passed:
            return self.llm.report("❌ QC failed. Check reports.")
        
        # 6. Ask for downstream confirmation
        if intent.includes_downstream:
            confirm = self.llm.confirm_with_user(
                "QC passed. Proceed with downstream analysis?"
            )
            if confirm:
                downstream = pipeline.run_downstream(qc_result.output)
                return downstream
```

---

## Configuration Structure

### Pipeline Paths Configuration

```yaml
# configs/pipeline_paths.yaml

pipelines:
  rna-seq:
    production: /data_3tb/shared/rna-seq-pipeline
    development: /home/ygkim/ngs-pipeline/rna-seq-pipeline
    conda_env: rna-seq-pipeline
    workflow: Snakefile
    output_base: /data_3tb/shared/output
    downstream:
      differential_expression:
        enabled: true
        path: /home/ygkim/ngs-pipeline/RNA-Seq_DE_GO_analysis
        auto_bridge: true
      gene_ontology:
        integrated: true  # Part of DE analysis
    
  atac-seq:
    production: /data_Xtb/atac-seq-pipeline
    conda_env: atac-seq-pipeline
    workflow: Snakefile
    output_base: /data_Xtb/output
    downstream:
      peak_calling:
        enabled: true
        tool: macs2
        params: "--broad"
      motif_analysis:
        enabled: true
        tool: homer
      differential_accessibility:
        enabled: true
        tool: diffbind
    
  wgs:
    production: /data_Ytb/wgs-pipeline
    conda_env: wgs-pipeline
    workflow: Snakefile
    output_base: /data_Ytb/output
    downstream:
      variant_calling:
        enabled: true
        tool: gatk
        mode: germline
      annotation:
        enabled: true
        tool: annovar
        databases: [refGene, clinvar, gnomad]
      quality_control:
        coverage_threshold: 30
        variant_quality_threshold: 30
```

---

## Implementation Roadmap

### Phase 5-7 (Current): RNA-seq Agent
- ✅ Phase 4: E2E testing completed
- 🔄 Phase 5: Ollama installation
- ⏳ Phase 6: LLM agent integration test
- ⏳ Phase 7: Production deployment

### Phase 8: Data Detection (Q2 2026)
- Implement `DataTypeDetector`
- Create detection rules
- Test with real datasets
- Validate accuracy (>95% target)

### Phase 9: ATAC-seq Integration (Q3 2026)
- Develop ATAC-seq pipeline (based on RNA-seq template)
- Implement `ATACSeqPipeline` wrapper
- Test peak calling, motif analysis
- Integrate into unified agent

### Phase 10: WGS Integration (Q4 2026)
- Develop WGS pipeline
- Implement `WGSPipeline` wrapper
- Test variant calling, annotation
- Full multi-pipeline integration

### Phase 11: Unified Agent (2027)
- Deploy `/data_3tb/shared/ngs-agent/`
- Natural language interface for all types
- Production monitoring and logging
- User training and documentation

---

## Benefits of Phased Approach

### Current (RNA-seq only)
✅ **Low risk**: Single pipeline, well-tested  
✅ **Fast validation**: Immediate feedback  
✅ **Template creation**: Standard for other pipelines  
✅ **Learning**: User experience, LLM capabilities

### Future (Unified)
🚀 **Scalability**: Easy to add new data types  
🚀 **Consistency**: Same interface for all pipelines  
🚀 **Efficiency**: Shared resources (Ollama, tools)  
🚀 **Maintainability**: Centralized agent logic

---

## Technical Considerations

### Ollama Placement
- **Current**: System-wide installation (`/usr/local/bin/ollama`)
- **Models**: `/data_3tb/shared/ollama-models/` (large disk)
- **Service**: systemd, single instance serves all pipelines
- **API**: `localhost:11434` (local only, secure)

### Resource Management
- **Memory**: 12-16GB for llama3.1:8b
- **CPU**: 4-8 cores recommended
- **Storage**: ~5GB per model
- **Concurrency**: 1 model loaded at a time (for now)

### Security
- ✅ Input validation (current)
- ✅ Path whitelisting (current)
- 🔄 Pipeline isolation (future)
- 🔄 Audit logging (future)
- 🔄 Rate limiting (future)

---

## Success Metrics

### Phase 5-7 (RNA-seq)
- [ ] Ollama installed and stable
- [ ] Agent responds to natural language
- [ ] All tools callable via LLM
- [ ] Security validation working
- [ ] Production deployment documented

### Phase 8-11 (Unified)
- [ ] Data type detection >95% accuracy
- [ ] All 3 pipelines integrated
- [ ] Single command for any NGS data
- [ ] User satisfaction survey >4/5
- [ ] Processing time <2x manual workflow

---

## References

- **Current Implementation**: `scripts/standardization/llm_agent.py`
- **Security**: `docs/developer/SECURITY_GUIDE.md`
- **Ollama Setup**: `docs/developer/PHASE5_OLLAMA_SETUP.md`
- **Testing**: `docs/developer/PHASE4_E2E_TESTING.md`

---

**Note**: This is a living document. Update as we progress through phases and learn from RNA-seq agent implementation.
