# RNA-seq Pipeline QC Report 가이드

## 📊 개요

RNA-seq 파이프라인은 모든 분석이 완료된 후 자동으로 **HTML 형식의 QC (Quality Control) 리포트**를 생성합니다. 이 리포트는 파이프라인의 각 단계별 품질 지표를 시각화하여 연구자가 실험 결과의 신뢰성을 쉽게 평가할 수 있도록 돕습니다.

### 리포트 위치
```
results/qc_report.html
```

### 리포트 열람 방법
- 웹 브라우저에서 직접 열기 (Chrome, Firefox, Edge 등)
- 별도의 소프트웨어 설치 불필요
- 인터넷 연결 없이도 열람 가능 (standalone HTML)

---

## 🎯 리포트 구성 요소

QC 리포트는 다음 6가지 주요 섹션으로 구성됩니다:

### 1. **파이프라인 개요 (Pipeline Overview)**
- 분석 완료 날짜 및 시간
- 전체 샘플 수
- 파이프라인 버전 정보

### 2. **Adapter Trimming 통계 (Cutadapt Results)**

각 샘플별로 어댑터 제거 과정의 결과를 표시합니다.

**주요 지표:**
- **Total Read Pairs**: 입력된 전체 read pair 수
- **Passed Filters**: 품질 필터를 통과한 read pair 수
- **Too Short**: 너무 짧아서 제거된 read 수
- **Pass Rate (%)**: 통과율 (높을수록 좋음, 일반적으로 >90%)

**해석:**
- ✅ **정상**: Pass rate 90-98%
- ⚠️ **주의**: Pass rate 80-90% (어댑터 오염 또는 품질 저하 가능성)
- ❌ **문제**: Pass rate <80% (실험 문제 의심)

**📚 참고 문헌:**
- Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. *EMBnet.journal*, 17(1), 10-12. [DOI: 10.14806/ej.17.1.200](https://doi.org/10.14806/ej.17.1.200)
- Chen, C., et al. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*, 34(17), i884-i890. [DOI: 10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560)

---

### 3. **Alignment 품질 (STAR Alignment Quality)**

각 샘플이 reference genome에 얼마나 잘 정렬되었는지 보여줍니다.

**주요 지표:**
- **Input Reads**: 입력된 총 read 수
- **Uniquely Mapped**: 게놈의 한 곳에만 매핑된 read 수
- **Uniquely Mapped %**: 고유 매핑 비율 (가장 중요한 지표)
- **Multi-mapped**: 여러 곳에 매핑된 read 수
- **Multi-mapped %**: 다중 매핑 비율
- **Unmapped (too short)**: 매핑되지 않은 read 수

**시각화:**
- 각 샘플별 매핑률을 **진행 바(Progress Bar)**로 표시
- 색상 코딩:
  - 🟢 **초록색** (Uniquely Mapped): 신뢰도 높음
  - 🟡 **노란색** (Multi-mapped): 주의 필요
  - 🔴 **빨간색** (Unmapped): 매핑 실패

**해석:**
- ✅ **정상 (고품질)**: Uniquely Mapped > 85%
- ✅ **정상 (양호)**: Uniquely Mapped 70-85%
- ⚠️ **주의**: Uniquely Mapped 50-70% (샘플 품질 또는 reference genome 문제)
- ❌ **문제**: Uniquely Mapped < 50% (심각한 문제, 재실험 고려)

**샘플별 권장 기준:**

| 실험 유형 | 최소 Uniquely Mapped % | 이상적 범위 |
|-----------|------------------------|-------------|
| mRNA-seq (poly-A) | 70% | 80-95% |
| Total RNA-seq | 60% | 70-85% |
| Degraded RNA | 50% | 60-75% |

**📚 참고 문헌:**
- Dobin, A., et al. (2013). STAR: ultrafast universal RNA-seq aligner. *Bioinformatics*, 29(1), 15-21. [DOI: 10.1093/bioinformatics/bts635](https://doi.org/10.1093/bioinformatics/bts635)
- ENCODE Project Consortium (2012). An integrated encyclopedia of DNA elements in the human genome. *Nature*, 489(7414), 57-74. [DOI: 10.1038/nature11247](https://doi.org/10.1038/nature11247)
- Conesa, A., et al. (2016). A survey of best practices for RNA-seq data analysis. *Genome Biology*, 17(1), 13. [DOI: 10.1186/s13059-016-0881-8](https://doi.org/10.1186/s13059-016-0881-8)

---

### 4. **Gene Quantification 통계 (featureCounts Summary)**

featureCounts 단계에서 각 샘플의 read가 유전자에 할당된 결과를 보여줍니다.

**주요 지표:**
- **Assigned**: 유전자에 성공적으로 할당된 read 수 (가장 중요)
- **Unassigned_MultiMapping**: 여러 유전자에 매핑되어 할당되지 않은 read
- **Unassigned_NoFeatures**: 유전자가 아닌 영역에 매핑된 read (intergenic/intronic)
- **Unassigned_Ambiguity**: 여러 유전자에 걸쳐 애매하게 매핑된 read

**해석:**
- ✅ **정상**: Assigned 60-80% (mRNA-seq)
- ✅ **정상**: Assigned 40-60% (Total RNA-seq)
- ⚠️ **주의**: Unassigned_NoFeatures > 30% (게놈 오염 또는 rRNA 제거 실패)

**Assignment Rate 계산:**
```
Assignment Rate (%) = (Assigned / Total Reads) × 100
```

**📚 참고 문헌:**
- Liao, Y., Smyth, G. K., & Shi, W. (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics*, 30(7), 923-930. [DOI: 10.1093/bioinformatics/btt656](https://doi.org/10.1093/bioinformatics/btt656)
- Liao, Y., Smyth, G. K., & Shi, W. (2019). The R package Rsubread is easier, faster, cheaper and better for alignment and quantification of RNA sequencing reads. *Nucleic Acids Research*, 47(8), e47. [DOI: 10.1093/nar/gkz114](https://doi.org/10.1093/nar/gkz114)

---

### 5. **유전자 발현 통계 (Gene Expression Statistics)**

각 샘플에서 검출된 유전자와 총 발현량에 대한 통계를 제공합니다.

**주요 지표:**
- **Total Counts**: 샘플의 총 read count 수
- **Genes Detected**: 발현이 검출된 유전자 수 (count > 0)

**샘플별 비교:**
- 샘플 간 Total Counts가 크게 차이나면 normalization 필요
- Genes Detected는 샘플 품질의 지표 (일반적으로 15,000-20,000개)

**해석 (Mouse 기준):**
- ✅ **정상**: Genes Detected 15,000-20,000개
- ⚠️ **주의**: Genes Detected 10,000-15,000개 (낮은 깊이 또는 품질)
- ❌ **문제**: Genes Detected < 10,000개 (심각한 문제)

**해석 (Human 기준):**
- ✅ **정상**: Genes Detected 15,000-25,000개
- ⚠️ **주의**: Genes Detected 12,000-15,000개
- ❌ **문제**: Genes Detected < 12,000개

**📚 참고 문헌:**
- Wang, Z., Gerstein, M., & Snyder, M. (2009). RNA-Seq: a revolutionary tool for transcriptomics. *Nature Reviews Genetics*, 10(1), 57-63. [DOI: 10.1038/nrg2484](https://doi.org/10.1038/nrg2484)
- Griffith, M., et al. (2015). Informatics for RNA sequencing: a web resource for analysis on the cloud. *PLoS Computational Biology*, 11(8), e1004393. [DOI: 10.1371/journal.pcbi.1004393](https://doi.org/10.1371/journal.pcbi.1004393)

---

### 6. **고발현 유전자 Top N (Top Expressed Genes)**

각 샘플에서 가장 많이 발현된 유전자 목록입니다.

**표시 내용:**
- Gene ID (Ensembl ID)
- Raw count 값
- 샘플별 상위 N개 유전자 (기본값: 10개)

**확인 사항:**
- **정상적인 경우**: 조직/세포 특이적 마커 유전자가 상위권에 있어야 함
- **문제 징후**:
  - Mitochondrial genes (MT-*) 과다 → 세포 스트레스/사멸
  - Ribosomal genes (RPL*, RPS*) 과다 → rRNA 제거 실패
  - Hemoglobin genes (HBA*, HBB*) 과다 → 혈액 오염

**일반적인 housekeeping genes:**
- ACTB, GAPDH, B2M, PPIA, RPLP0 등이 중간-높은 발현 보임

**📚 참고 문헌:**
- Eisenberg, E., & Levanon, E. Y. (2013). Human housekeeping genes, revisited. *Trends in Genetics*, 29(10), 569-574. [DOI: 10.1016/j.tig.2013.05.010](https://doi.org/10.1016/j.tig.2013.05.010)
- de Jonge, H. J., et al. (2007). Evidence based selection of housekeeping genes. *PloS One*, 2(9), e898. [DOI: 10.1371/journal.pone.0000898](https://doi.org/10.1371/journal.pone.0000898)
- Sheng, Q., et al. (2017). Multi-perspective quality control of Illumina RNA sequencing data analysis. *Briefings in Functional Genomics*, 16(4), 194-204. [DOI: 10.1093/bfgp/elw035](https://doi.org/10.1093/bfgp/elw035)

---

## ⚙️ 리포트 설정 및 커스터마이징

### 1. QC 리포트 생성 활성화/비활성화

`config.yaml` 파일에서 설정:

```yaml
# === QC Report Parameters ===
generate_qc_report: true   # true: 리포트 생성, false: 비활성화
qc_report_output: "results/qc_report.html"  # 출력 파일 경로
qc_top_genes: 10  # 상위 발현 유전자 표시 개수
```

### 2. 조정 가능한 파라미터

| 파라미터 | 설명 | 기본값 | 권장 범위 |
|---------|------|--------|----------|
| `generate_qc_report` | 리포트 생성 여부 | `true` | true/false |
| `qc_report_output` | 출력 파일 경로 | `"results/qc_report.html"` | 원하는 경로 |
| `qc_top_genes` | 상위 유전자 표시 개수 | `10` | 5-20 |

### 3. Top Genes 개수 변경 예시

**더 많은 유전자 확인 (20개):**
```yaml
qc_top_genes: 20
```

**적은 개수로 간단히 (5개):**
```yaml
qc_top_genes: 5
```

### 4. 리포트 재생성

설정을 변경한 후 리포트만 다시 생성하려면:

```bash
cd /home/ygkim/ngs_pipeline/rna-seq-pipeline

# 리포트만 강제 재생성
snakemake --cores 1 results/qc_report.html --force
```

---

## 📈 QC 리포트 활용 방법

### 1. **실험 품질 평가**
- 모든 샘플이 유사한 품질 지표를 보이는지 확인
- 이상치(outlier) 샘플 식별

### 2. **샘플 간 비교**
- Total counts, mapping rate, genes detected 비교
- 배치 효과(batch effect) 확인

### 3. **문제 진단**
- 낮은 mapping rate → reference genome 문제 또는 오염
- 낮은 assignment rate → annotation 파일 문제
- 적은 genes detected → 시퀀싱 깊이 부족 또는 RNA 품질 저하

### 4. **다운스트림 분석 판단**
- 품질이 좋은 샘플 선별
- 문제 샘플 제외 결정
- 추가 시퀀싱 필요성 판단

**📚 참고 문헌:**
- García-Alcalde, F., et al. (2012). Qualimap: evaluating next-generation sequencing alignment data. *Bioinformatics*, 28(20), 2678-2679. [DOI: 10.1093/bioinformatics/bts503](https://doi.org/10.1093/bioinformatics/bts503)
- Ewels, P., et al. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics*, 32(19), 3047-3048. [DOI: 10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354)

---

## 🔍 샘플 품질 체크리스트

각 샘플이 다음 기준을 충족하는지 확인하세요:

### ✅ 우수한 품질
- [ ] Cutadapt pass rate > 90%
- [ ] Uniquely mapped % > 85%
- [ ] Assignment rate > 70%
- [ ] Genes detected > 15,000개

### ⚠️ 양호한 품질 (사용 가능)
- [ ] Cutadapt pass rate > 85%
- [ ] Uniquely mapped % > 70%
- [ ] Assignment rate > 60%
- [ ] Genes detected > 12,000개

### ❌ 문제가 있는 품질 (재검토 필요)
- [ ] Cutadapt pass rate < 80%
- [ ] Uniquely mapped % < 60%
- [ ] Assignment rate < 50%
- [ ] Genes detected < 10,000개

---

## 📋 리포트 공유 및 보관

### 1. 보고서 첨부
- HTML 파일을 이메일이나 공유 폴더에 첨부
- 브라우저만 있으면 누구나 열람 가능

### 2. 스크린샷 활용
- 주요 테이블이나 차트를 캡처하여 논문/보고서에 삽입
- 브라우저의 개발자 도구로 특정 섹션만 출력 가능

### 3. 데이터 보관
- 분석 결과와 함께 QC 리포트 보관 (재현성 확보)
- 파일 경로: `results/qc_report.html`

---

## 🛠️ 고급 활용

### 1. 스크립트 직접 실행

Snakemake 외부에서 직접 QC 리포트 생성:

```bash
cd /home/ygkim/ngs_pipeline/rna-seq-pipeline
python src/generate_qc_report.py
```

### 2. 커스텀 분석 추가

`src/generate_qc_report.py` 스크립트를 수정하여 추가 분석 가능:
- 샘플 간 상관관계 분석
- PCA plot 추가
- 커스텀 품질 지표 계산

### 3. 리포트 스타일 변경

HTML/CSS를 수정하여 디자인 커스터마이징 가능:
- 회사/연구실 로고 추가
- 색상 테마 변경
- 추가 설명 섹션 삽입

---

## 📞 문제 해결

### Q1. QC 리포트가 생성되지 않아요
**A:** `config.yaml`에서 `generate_qc_report: true`로 설정되어 있는지 확인하세요.

### Q2. 리포트에 일부 샘플 데이터가 누락되었어요
**A:** 해당 샘플의 로그 파일이나 결과 파일이 제대로 생성되었는지 확인:
- `logs/cutadapt/{sample}.log`
- `results/aligned/{sample}/Log.final.out`
- `results/counts/counts_matrix.txt`

### Q3. Top genes 목록에 Gene Symbol이 아닌 ID만 나와요
**A:** 현재 버전은 Ensembl ID만 표시합니다. Gene symbol 변환 기능은 추후 업데이트 예정입니다.

### Q4. 리포트를 PDF로 변환하고 싶어요
**A:** 브라우저의 인쇄 기능(Ctrl+P 또는 Cmd+P)을 사용하여 PDF로 저장할 수 있습니다.

---

## 📚 종합 참고 문헌

### RNA-seq 품질 관리 종합 리뷰
1. **Conesa, A., et al. (2016).** A survey of best practices for RNA-seq data analysis. *Genome Biology*, 17(1), 13. [DOI: 10.1186/s13059-016-0881-8](https://doi.org/10.1186/s13059-016-0881-8)
   - RNA-seq 분석의 모든 단계에 대한 종합적인 best practice 가이드

2. **Williams, C. R., et al. (2017).** Empirical assessment of analysis workflows for differential expression analysis of human samples using RNA-Seq. *BMC Bioinformatics*, 18(1), 38. [DOI: 10.1186/s12859-016-1457-z](https://doi.org/10.1186/s12859-016-1457-z)
   - 다양한 RNA-seq 분석 워크플로우의 비교 평가

3. **Costa-Silva, J., Domingues, D., & Lopes, F. M. (2017).** RNA-Seq differential expression analysis: An extended review and a software tool. *PloS One*, 12(12), e0190152. [DOI: 10.1371/journal.pone.0190152](https://doi.org/10.1371/journal.pone.0190152)
   - RNA-seq 차등 발현 분석의 확장된 리뷰

### 개별 도구 원문
4. **Martin, M. (2011).** Cutadapt removes adapter sequences from high-throughput sequencing reads. *EMBnet.journal*, 17(1), 10-12. [DOI: 10.14806/ej.17.1.200](https://doi.org/10.14806/ej.17.1.200)

5. **Dobin, A., et al. (2013).** STAR: ultrafast universal RNA-seq aligner. *Bioinformatics*, 29(1), 15-21. [DOI: 10.1093/bioinformatics/bts635](https://doi.org/10.1093/bioinformatics/bts635)

6. **Liao, Y., Smyth, G. K., & Shi, W. (2014).** featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics*, 30(7), 923-930. [DOI: 10.1093/bioinformatics/btt656](https://doi.org/10.1093/bioinformatics/btt656)

### 품질 평가 기준
7. **ENCODE Project Consortium (2012).** An integrated encyclopedia of DNA elements in the human genome. *Nature*, 489(7414), 57-74. [DOI: 10.1038/nature11247](https://doi.org/10.1038/nature11247)
   - ENCODE 프로젝트의 RNA-seq 품질 기준

8. **Wang, L., et al. (2016).** Measure transcript integrity using RNA-seq data. *BMC Bioinformatics*, 17(1), 58. [DOI: 10.1186/s12859-016-0922-z](https://doi.org/10.1186/s12859-016-0922-z)
   - RNA 품질 평가 방법론

9. **Sheng, Q., et al. (2017).** Multi-perspective quality control of Illumina RNA sequencing data analysis. *Briefings in Functional Genomics*, 16(4), 194-204. [DOI: 10.1093/bfgp/elw035](https://doi.org/10.1093/bfgp/elw035)
   - 다각적 RNA-seq QC 방법론

### Housekeeping Genes 및 바이오마커
10. **Eisenberg, E., & Levanon, E. Y. (2013).** Human housekeeping genes, revisited. *Trends in Genetics*, 29(10), 569-574. [DOI: 10.1016/j.tig.2013.05.010](https://doi.org/10.1016/j.tig.2013.05.010)

11. **de Jonge, H. J., et al. (2007).** Evidence based selection of housekeeping genes. *PloS One*, 2(9), e898. [DOI: 10.1371/journal.pone.0000898](https://doi.org/10.1371/journal.pone.0000898)

### 온라인 자료 및 가이드라인
- **ENCODE RNA-seq Standards**: https://www.encodeproject.org/data-standards/rna-seq/
- **RNA-seqlopedia**: https://rnaseq.uoregon.edu/ (오레곤 대학교 RNA-seq 교육 자료)
- **Cutadapt Documentation**: https://cutadapt.readthedocs.io/
- **STAR Manual**: https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
- **Subread/featureCounts**: http://subread.sourceforge.net/
- **RNAseq Analysis Guide (Harvard Chan Bioinformatics Core)**: https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon/

### 추가 QC 도구 참고 문헌
12. **Andrews, S. (2010).** FastQC: a quality control tool for high throughput sequence data. Available online at: http://www.bioinformatics.babraham.ac.uk/projects/fastqc

13. **Ewels, P., et al. (2016).** MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics*, 32(19), 3047-3048. [DOI: 10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354)

14. **Wang, L., et al. (2012).** RSeQC: quality control of RNA-seq experiments. *Bioinformatics*, 28(16), 2184-2185. [DOI: 10.1093/bioinformatics/bts356](https://doi.org/10.1093/bioinformatics/bts356)

---

## 📝 업데이트 이력

- **2025-12-01**: 초기 문서 작성
  - QC 리포트 구성 요소 상세 설명
  - 품질 평가 기준 추가
  - 커스터마이징 가이드 추가

---

**문의 사항이나 개선 제안이 있으시면 이슈를 등록해주세요!** 🙌
