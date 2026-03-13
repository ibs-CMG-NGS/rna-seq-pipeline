#!/usr/bin/env python3
"""
RNA-seq Pipeline QC Report Generator
여러 샘플의 결과를 종합하여 HTML 리포트 생성
Snakemake script로 사용 가능
"""

import os
import glob
import json
from datetime import datetime
from pathlib import Path

# Snakemake에서 실행될 때 사용
try:
    # snakemake 객체가 있으면 script mode
    TOP_GENES = snakemake.params.get('top_genes', 10)
    OUTPUT_FILE = str(snakemake.output[0])
    LOG_FILE = str(snakemake.log[0])
    SCRIPT_MODE = True
    
    # Snakemake params에서 절대 경로 직접 가져오기
    LOGS_DIR = snakemake.params.get('logs_dir')
    RESULTS_DIR = snakemake.params.get('results_dir')
    DATA_DIR = snakemake.params.get('data_dir')
    BASE_DIR = DATA_DIR  # data_dir이 base directory
    TRIMMED_DIR = snakemake.params.get('trimmed_dir', os.path.join(RESULTS_DIR, 'trimmed'))
    ALIGNED_DIR = snakemake.params.get('aligned_dir', os.path.join(RESULTS_DIR, 'aligned'))
    COUNTS_DIR  = snakemake.params.get('counts_dir', os.path.join(RESULTS_DIR, 'counts'))
    RAW_FILES   = snakemake.params.get('raw_files', {})
    PIPELINE_QC_JSON = str(snakemake.input.pipeline_qc_json)

except NameError:
    # 직접 실행 모드
    TOP_GENES = 10
    OUTPUT_FILE = 'results/qc_report.html'
    LOG_FILE = None
    SCRIPT_MODE = False
    RESULTS_DIR = 'results'
    BASE_DIR = '.'
    LOGS_DIR = 'logs'
    DATA_DIR = 'data'
    TRIMMED_DIR = os.path.join(RESULTS_DIR, 'trimmed')
    ALIGNED_DIR = os.path.join(RESULTS_DIR, 'aligned')
    COUNTS_DIR  = os.path.join(RESULTS_DIR, 'counts')
    RAW_FILES   = {}
    PIPELINE_QC_JSON = None

def get_sample_names():
    """trimmed 폴더에서 샘플 이름 추출"""
    samples = set()
    pattern = os.path.join(TRIMMED_DIR, '*_1.fastq.gz')
    for f in glob.glob(pattern):
        sample = os.path.basename(f).replace('_1.fastq.gz', '')
        samples.add(sample)
    if not samples:
        print(f"❌ No samples found in {TRIMMED_DIR}")
    return sorted(samples)

def parse_cutadapt_log(sample):
    """Cutadapt 로그 파싱"""
    log_file = os.path.join(LOGS_DIR, 'cutadapt', f'{sample}.log')
    if not os.path.exists(log_file):
        return None
    
    import re
    data = {}
    with open(log_file, 'r') as f:
        content = f.read()
        for line in content.split('\n'):
            if 'Total read pairs processed' in line:
                # Extract number with commas: "52,280,544"
                m = re.search(r'([\d,]+)', line)
                if m:
                    data['total_pairs'] = int(m.group(1).replace(',', ''))
            elif 'Pairs that were too short' in line:
                m = re.search(r'([\d,]+)', line)
                if m:
                    data['too_short'] = int(m.group(1).replace(',', ''))
            elif 'Pairs written (passing filters)' in line:
                m = re.search(r'([\d,]+)', line)
                if m:
                    data['passed'] = int(m.group(1).replace(',', ''))
            elif 'Total basepairs processed' in line:
                m = re.search(r'([\d,]+)', line)
                if m:
                    data['total_bp'] = int(m.group(1).replace(',', ''))
    return data

def parse_star_log(sample):
    """STAR Log.final.out 파싱"""
    log_file = os.path.join(ALIGNED_DIR, sample, 'Log.final.out')
    if not os.path.exists(log_file):
        return None
    
    import re
    data = {}
    with open(log_file, 'r') as f:
        for line in f:
            line = line.strip()
            if 'Number of input reads' in line:
                m = re.search(r'(\d[\d,]*)', line)
                if m:
                    data['input_reads'] = int(m.group(1).replace(',', ''))
            elif 'Uniquely mapped reads number' in line:
                m = re.search(r'(\d[\d,]*)', line)
                if m:
                    data['unique_mapped'] = int(m.group(1).replace(',', ''))
            elif 'Uniquely mapped reads %' in line:
                m = re.search(r'([\d.]+)%', line)
                if m:
                    data['unique_mapped_pct'] = float(m.group(1))
            elif 'Number of reads mapped to multiple loci' in line:
                m = re.search(r'(\d[\d,]*)', line)
                if m:
                    data['multi_mapped'] = int(m.group(1).replace(',', ''))
            elif '% of reads mapped to multiple loci' in line:
                m = re.search(r'([\d.]+)%', line)
                if m:
                    data['multi_mapped_pct'] = float(m.group(1))
            elif 'Number of reads unmapped: too short' in line:
                m = re.search(r'(\d[\d,]*)', line)
                if m:
                    data['unmapped_short'] = int(m.group(1).replace(',', ''))
            elif '% of reads unmapped: too short' in line:
                m = re.search(r'([\d.]+)%', line)
                if m:
                    data['unmapped_short_pct'] = float(m.group(1))
    return data

def parse_featurecounts_summary():
    """featureCounts summary 파싱"""
    summary_file = os.path.join(RESULTS_DIR, 'counts', 'counts_matrix.txt.summary')
    if not os.path.exists(summary_file):
        return None
    
    data = {}
    with open(summary_file, 'r') as f:
        lines = f.readlines()
        header = lines[0].strip().split('\t')
        
        # 샘플별 데이터
        for col_idx in range(1, len(header)):
            sample_path = header[col_idx]
            sample_name = sample_path.split('/')[-2]  # aligned/SAMPLE/file.bam
            data[sample_name] = {}
            
            for line in lines[1:]:
                parts = line.strip().split('\t')
                if len(parts) > col_idx:
                    status = parts[0]
                    count = int(parts[col_idx])
                    data[sample_name][status] = count
    
    return data

def analyze_count_matrix():
    """Count matrix 분석"""
    count_file = os.path.join(RESULTS_DIR, 'counts', 'counts_matrix.txt')
    if not os.path.exists(count_file):
        return None
    
    samples_data = {}
    gene_stats = {'total_genes': 0, 'genes_with_reads': 0}
    
    with open(count_file, 'r') as f:
        for i, line in enumerate(f):
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            
            # 헤더 라인 (샘플명 추출)
            if i == 1:
                sample_columns = {}
                for col_idx in range(6, len(parts)):
                    sample_path = parts[col_idx]
                    sample_name = sample_path.split('/')[-2]
                    sample_columns[col_idx] = sample_name
                    samples_data[sample_name] = {
                        'total_counts': 0,
                        'genes_detected': 0,
                        'high_expression': []
                    }
                continue
            
            # 데이터 라인
            gene_id = parts[0]
            gene_stats['total_genes'] += 1
            
            for col_idx, sample_name in sample_columns.items():
                count = int(parts[col_idx])
                samples_data[sample_name]['total_counts'] += count
                
                if count > 0:
                    samples_data[sample_name]['genes_detected'] += 1
                
                if count > 10000:
                    samples_data[sample_name]['high_expression'].append((gene_id, count))
    
    # 상위 발현 유전자 정렬
    for sample in samples_data:
        samples_data[sample]['high_expression'].sort(key=lambda x: x[1], reverse=True)
        samples_data[sample]['high_expression'] = samples_data[sample]['high_expression'][:TOP_GENES]
    
    return {'samples': samples_data, 'gene_stats': gene_stats}

def _render_pipeline_qc_section(json_path):
    """Pipeline QC evaluation JSON을 읽어 HTML 섹션으로 렌더링 (English)."""
    import json as _json

    if not json_path or not os.path.exists(json_path):
        return ""

    try:
        with open(json_path, encoding='utf-8') as f:
            data = _json.load(f)
    except Exception:
        return ""

    samples = data.get('samples', [])
    thresholds = data.get('thresholds', {})
    if not samples:
        return ""

    status_color = {'PASS': '#10b981', 'WARN': '#f59e0b', 'FAIL': '#ef4444'}
    status_bg    = {'PASS': '#d1fae5', 'WARN': '#fef3c7', 'FAIL': '#fee2e2'}
    status_icon  = {'PASS': '✅', 'WARN': '⚠️', 'FAIL': '❌'}

    counts = {'PASS': 0, 'WARN': 0, 'FAIL': 0}
    for s in samples:
        counts[s['status']] = counts.get(s['status'], 0) + 1

    html = """
            <!-- Pipeline QC Evaluation Summary -->
            <div class="section">
                <h2 class="section-title">🔬 Pipeline QC Evaluation Summary</h2>
                <p style="color:#555; margin-bottom:16px;">
                    Per-sample judgment based on Cutadapt, STAR, and featureCounts metrics.
                </p>
"""

    # Summary cards
    html += f"""
                <div class="summary-grid" style="grid-template-columns: repeat(3, 1fr); max-width:600px; margin-bottom:24px;">
                    <div class="summary-card" style="background:{status_bg['PASS']};">
                        <h3 style="color:{status_color['PASS']};">PASS</h3>
                        <div class="value" style="color:{status_color['PASS']};">{counts['PASS']}</div>
                        <div class="sub-value">Ready for DE analysis</div>
                    </div>
                    <div class="summary-card" style="background:{status_bg['WARN']};">
                        <h3 style="color:{status_color['WARN']};">WARN</h3>
                        <div class="value" style="color:{status_color['WARN']};">{counts['WARN']}</div>
                        <div class="sub-value">Proceed with caution</div>
                    </div>
                    <div class="summary-card" style="background:{status_bg['FAIL']};">
                        <h3 style="color:{status_color['FAIL']};">FAIL</h3>
                        <div class="value" style="color:{status_color['FAIL']};">{counts['FAIL']}</div>
                        <div class="sub-value">Review before analysis</div>
                    </div>
                </div>
"""

    # Thresholds note
    html += f"""
                <details style="margin-bottom:16px; background:#f9fafb; padding:12px; border-radius:8px;">
                    <summary style="cursor:pointer; font-weight:600; color:#667eea;">Applied Thresholds</summary>
                    <ul style="margin-top:8px; color:#555; font-size:0.9em;">
                        <li>Cutadapt too-short: WARN ≥ {thresholds.get('cutadapt_too_short_warn',2)}% / FAIL ≥ {thresholds.get('cutadapt_too_short_fail',5)}%</li>
                        <li>STAR uniquely mapped: WARN &lt; {thresholds.get('min_uniquely_mapped_pct_warn',75)}% / FAIL &lt; {thresholds.get('min_uniquely_mapped_pct_fail',70)}%</li>
                        <li>STAR unmapped (too short): WARN &gt; {thresholds.get('max_unmapped_short_warn',13)}% / FAIL &gt; {thresholds.get('max_unmapped_short_fail',17)}%</li>
                        <li>STAR mismatch rate: WARN &gt; {thresholds.get('max_mismatch_rate_warn',1.5)}% / FAIL &gt; {thresholds.get('max_mismatch_rate_fail',2)}%</li>
                        <li>featureCounts assigned: WARN &lt; {thresholds.get('min_assignment_rate_warn',50)}% / FAIL &lt; {thresholds.get('min_assignment_rate_fail',30)}%</li>
                    </ul>
                </details>
"""

    # Per-sample table
    html += """
                <table>
                    <thead>
                        <tr>
                            <th>Sample</th>
                            <th>Status</th>
                            <th>Cutadapt<br>too-short</th>
                            <th>STAR<br>unique %</th>
                            <th>STAR<br>unmapped %</th>
                            <th>STAR<br>mismatch %</th>
                            <th>featureCounts<br>assigned %</th>
                            <th>Issues / Warnings</th>
                        </tr>
                    </thead>
                    <tbody>
"""
    for s in sorted(samples, key=lambda x: (x['status'] != 'FAIL', x['status'] != 'WARN', x['sample'])):
        m = s.get('metrics', {})
        bg = status_bg[s['status']]
        color = status_color[s['status']]
        icon = status_icon[s['status']]

        def fmt(key, fmt_str='{:.1f}%'):
            v = m.get(key)
            return fmt_str.format(v) if v is not None else 'N/A'

        notes = []
        for iss in s.get('issues', []):
            notes.append(f'<span style="color:{status_color["FAIL"]}; font-size:0.85em;">● {iss}</span>')
        for w in s.get('warnings', []):
            notes.append(f'<span style="color:{status_color["WARN"]}; font-size:0.85em;">○ {w}</span>')
        notes_html = '<br>'.join(notes) if notes else '—'

        html += f"""
                        <tr style="background:{bg}20;">
                            <td><strong>{s['sample']}</strong></td>
                            <td style="text-align:center;">
                                <span style="background:{bg}; color:{color}; padding:3px 10px; border-radius:10px; font-weight:600; font-size:0.85em;">
                                    {icon} {s['status']}
                                </span>
                            </td>
                            <td style="text-align:right;">{fmt('cutadapt_too_short_pct')}</td>
                            <td style="text-align:right;">{fmt('star_unique_mapped_pct', '{:.2f}%')}</td>
                            <td style="text-align:right;">{fmt('star_unmapped_short_pct', '{:.2f}%')}</td>
                            <td style="text-align:right;">{fmt('star_mismatch_rate', '{:.2f}%')}</td>
                            <td style="text-align:right;">{fmt('featurecounts_assignment_pct')}</td>
                            <td style="font-size:0.85em; line-height:1.6;">{notes_html}</td>
                        </tr>
"""

    html += """
                    </tbody>
                </table>
            </div>
"""
    return html


def get_file_size(filepath):
    """파일 크기를 읽기 쉬운 형식으로 변환"""
    if not os.path.exists(filepath):
        return "N/A"
    size = os.path.getsize(filepath)
    if size > 1e9:
        return f"{size/1e9:.2f} GB"
    else:
        return f"{size/1e6:.2f} MB"

def generate_html_report(samples):
    """HTML 리포트 생성"""
    
    # 각 샘플별 데이터 수집
    all_cutadapt = {}
    all_star = {}
    
    for sample in samples:
        all_cutadapt[sample] = parse_cutadapt_log(sample)
        all_star[sample] = parse_star_log(sample)
    
    fc_summary = parse_featurecounts_summary()
    count_analysis = analyze_count_matrix()
    
    # HTML 생성
    html = f"""
<!DOCTYPE html>
<html lang="ko">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>RNA-seq Pipeline QC Report</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            padding: 20px;
            color: #333;
        }}
        
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            background: white;
            border-radius: 15px;
            box-shadow: 0 10px 40px rgba(0,0,0,0.2);
            overflow: hidden;
        }}
        
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 40px;
            text-align: center;
        }}
        
        .header h1 {{
            font-size: 2.5em;
            margin-bottom: 10px;
        }}
        
        .header p {{
            font-size: 1.1em;
            opacity: 0.9;
        }}
        
        .content {{
            padding: 40px;
        }}
        
        .section {{
            margin-bottom: 40px;
        }}
        
        .section-title {{
            font-size: 1.8em;
            color: #667eea;
            margin-bottom: 20px;
            padding-bottom: 10px;
            border-bottom: 3px solid #667eea;
        }}
        
        .summary-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}
        
        .summary-card {{
            background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
            padding: 25px;
            border-radius: 10px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }}
        
        .summary-card h3 {{
            color: #667eea;
            font-size: 0.9em;
            margin-bottom: 10px;
            text-transform: uppercase;
        }}
        
        .summary-card .value {{
            font-size: 2em;
            font-weight: bold;
            color: #333;
        }}
        
        .summary-card .sub-value {{
            font-size: 0.9em;
            color: #666;
            margin-top: 5px;
        }}
        
        table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 20px;
            background: white;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            border-radius: 8px;
            overflow: hidden;
        }}
        
        th {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 15px;
            text-align: left;
            font-weight: 600;
        }}
        
        td {{
            padding: 12px 15px;
            border-bottom: 1px solid #eee;
        }}
        
        tr:hover {{
            background-color: #f8f9ff;
        }}
        
        .metric-good {{
            color: #10b981;
            font-weight: bold;
        }}
        
        .metric-warning {{
            color: #f59e0b;
            font-weight: bold;
        }}
        
        .metric-bad {{
            color: #ef4444;
            font-weight: bold;
        }}
        
        .progress-bar {{
            height: 25px;
            background: #e5e7eb;
            border-radius: 12px;
            overflow: hidden;
            margin: 5px 0;
        }}
        
        .progress-fill {{
            height: 100%;
            background: linear-gradient(90deg, #667eea 0%, #764ba2 100%);
            display: flex;
            align-items: center;
            justify-content: center;
            color: white;
            font-size: 0.85em;
            font-weight: bold;
            transition: width 0.3s ease;
        }}
        
        .gene-list {{
            background: #f9fafb;
            padding: 15px;
            border-radius: 8px;
            margin-top: 10px;
        }}
        
        .gene-item {{
            display: flex;
            justify-content: space-between;
            padding: 8px 0;
            border-bottom: 1px solid #e5e7eb;
        }}
        
        .gene-item:last-child {{
            border-bottom: none;
        }}
        
        .footer {{
            background: #f9fafb;
            padding: 20px;
            text-align: center;
            color: #666;
            font-size: 0.9em;
        }}
        
        .badge {{
            display: inline-block;
            padding: 4px 12px;
            border-radius: 12px;
            font-size: 0.85em;
            font-weight: 600;
        }}
        
        .badge-success {{
            background: #d1fae5;
            color: #065f46;
        }}
        
        .badge-warning {{
            background: #fef3c7;
            color: #92400e;
        }}
        
        .badge-info {{
            background: #dbeafe;
            color: #1e40af;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🧬 RNA-seq Pipeline QC Report</h1>
            <p>Quality Control & Analysis Summary</p>
            <p style="font-size: 0.9em; margin-top: 10px;">Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        </div>
        
        <div class="content">
            <!-- 전체 요약 -->
            <div class="section">
                <h2 class="section-title">📊 Overall Summary</h2>
                <div class="summary-grid">
                    <div class="summary-card">
                        <h3>Total Samples</h3>
                        <div class="value">{len(samples)}</div>
                        <div class="sub-value">{', '.join(samples)}</div>
                    </div>
"""
    
    # 전체 통계 추가
    if count_analysis:
        total_genes = count_analysis['gene_stats']['total_genes']
        total_counts = sum(s['total_counts'] for s in count_analysis['samples'].values())
        avg_genes_detected = sum(s['genes_detected'] for s in count_analysis['samples'].values()) / len(samples)
        
        html += f"""
                    <div class="summary-card">
                        <h3>Total Genes</h3>
                        <div class="value">{total_genes:,}</div>
                        <div class="sub-value">Avg detected: {avg_genes_detected:,.0f}</div>
                    </div>
                    <div class="summary-card">
                        <h3>Total Counts</h3>
                        <div class="value">{total_counts:,}</div>
                        <div class="sub-value">All samples combined</div>
                    </div>
"""
    
    html += """
                </div>
            </div>
            
            <!-- Cutadapt 결과 -->
            <div class="section">
                <h2 class="section-title">✂️ Adapter Trimming (Cutadapt)</h2>
                <table>
                    <thead>
                        <tr>
                            <th>Sample</th>
                            <th>Total Pairs</th>
                            <th>Passed</th>
                            <th>Too Short</th>
                            <th>Pass Rate</th>
                        </tr>
                    </thead>
                    <tbody>
"""
    
    for sample in samples:
        data = all_cutadapt.get(sample)
        if data:
            html += f"""
                        <tr>
                            <td><strong>{sample}</strong></td>
                            <td>{data.get('total_pairs', 'N/A'):,}</td>
                            <td>{data.get('passed', 'N/A'):,}</td>
                            <td>{data.get('too_short', 'N/A'):,}</td>
                            <td>
"""
            # Pass rate 계산 및 진행 바 추가
            if 'passed' in data and 'total_pairs' in data:
                try:
                    passed = data['passed']
                    total = data['total_pairs']
                    pass_rate = (passed / total) * 100
                    color_class = 'metric-good' if pass_rate > 95 else 'metric-warning'
                    html += f"""
                                <span class="{color_class}">{pass_rate:.1f}%</span>
                                <div class="progress-bar">
                                    <div class="progress-fill" style="width: {pass_rate}%">{pass_rate:.1f}%</div>
                                </div>
"""
                except:
                    html += "N/A"
            else:
                html += "N/A"
            
            html += """
                            </td>
                        </tr>
"""
    
    html += """
                    </tbody>
                </table>
            </div>
            
            <!-- STAR 정렬 결과 -->
            <div class="section">
                <h2 class="section-title">🎯 Alignment (STAR)</h2>
                <table>
                    <thead>
                        <tr>
                            <th>Sample</th>
                            <th>Input Reads</th>
                            <th>Uniquely Mapped</th>
                            <th>Multi-mapped</th>
                            <th>Unmapped (short)</th>
                            <th>Mapping Rate</th>
                        </tr>
                    </thead>
                    <tbody>
"""
    
    for sample in samples:
        data = all_star.get(sample)
        if data:
            unique_pct = data.get('unique_mapped_pct', 'N/A')
            multi_pct = data.get('multi_mapped_pct', 'N/A')
            unmapped_pct = data.get('unmapped_short_pct', 'N/A')
            
            # 매핑률 평가
            try:
                rate = float(unique_pct) if isinstance(unique_pct, (int, float)) else float(unique_pct.replace('%', ''))
                if rate > 85:
                    badge = '<span class="badge badge-success">Excellent</span>'
                    color_class = 'metric-good'
                elif rate > 70:
                    badge = '<span class="badge badge-warning">Good</span>'
                    color_class = 'metric-warning'
                else:
                    badge = '<span class="badge badge-warning">Check</span>'
                    color_class = 'metric-warning'
            except:
                badge = ''
                color_class = ''
            
            # Format numbers with commas
            input_reads = f"{data.get('input_reads', 'N/A'):,}" if isinstance(data.get('input_reads'), int) else data.get('input_reads', 'N/A')
            unique_mapped = f"{data.get('unique_mapped', 'N/A'):,}" if isinstance(data.get('unique_mapped'), int) else data.get('unique_mapped', 'N/A')
            multi_mapped = f"{data.get('multi_mapped', 'N/A'):,}" if isinstance(data.get('multi_mapped'), int) else data.get('multi_mapped', 'N/A')
            unmapped_short = f"{data.get('unmapped_short', 'N/A'):,}" if isinstance(data.get('unmapped_short'), int) else data.get('unmapped_short', 'N/A')
            
            # Format percentages
            unique_pct_str = f"({unique_pct:.2f}%)" if isinstance(unique_pct, (int, float)) else f"({unique_pct})"
            multi_pct_str = f"({multi_pct:.2f}%)" if isinstance(multi_pct, (int, float)) else f"({multi_pct})"
            unmapped_pct_str = f"({unmapped_pct:.2f}%)" if isinstance(unmapped_pct, (int, float)) else f"({unmapped_pct})"
            
            html += f"""
                        <tr>
                            <td><strong>{sample}</strong></td>
                            <td>{input_reads}</td>
                            <td>{unique_mapped} <span class="{color_class}">{unique_pct_str}</span></td>
                            <td>{multi_mapped} {multi_pct_str}</td>
                            <td>{unmapped_short} {unmapped_pct_str}</td>
                            <td>{badge}</td>
                        </tr>
"""
    
    html += """
                    </tbody>
                </table>
            </div>
            
            <!-- featureCounts 결과 -->
            <div class="section">
                <h2 class="section-title">🧮 Gene Quantification (featureCounts)</h2>
                <table>
                    <thead>
                        <tr>
                            <th>Sample</th>
                            <th>Assigned</th>
                            <th>MultiMapping</th>
                            <th>NoFeatures</th>
                            <th>Ambiguity</th>
                            <th>Assignment Rate</th>
                        </tr>
                    </thead>
                    <tbody>
"""
    
    if fc_summary:
        for sample in samples:
            data = fc_summary.get(sample, {})
            assigned = data.get('Assigned', 0)
            multi = data.get('Unassigned_MultiMapping', 0)
            no_features = data.get('Unassigned_NoFeatures', 0)
            ambiguity = data.get('Unassigned_Ambiguity', 0)
            total = sum(data.values())
            
            if total > 0:
                assign_rate = (assigned / total) * 100
                color_class = 'metric-good' if assign_rate > 70 else 'metric-warning'
            else:
                assign_rate = 0
                color_class = ''
            
            html += f"""
                        <tr>
                            <td><strong>{sample}</strong></td>
                            <td class="metric-good">{assigned:,}</td>
                            <td>{multi:,}</td>
                            <td>{no_features:,}</td>
                            <td>{ambiguity:,}</td>
                            <td>
                                <span class="{color_class}">{assign_rate:.1f}%</span>
                                <div class="progress-bar">
                                    <div class="progress-fill" style="width: {assign_rate}%">{assign_rate:.1f}%</div>
                                </div>
                            </td>
                        </tr>
"""
    
    html += """
                    </tbody>
                </table>
            </div>
            
            <!-- 유전자 발현 통계 -->
            <div class="section">
                <h2 class="section-title">📈 Gene Expression Statistics</h2>
                <table>
                    <thead>
                        <tr>
                            <th>Sample</th>
                            <th>Total Counts</th>
                            <th>Genes Detected</th>
                            <th>Detection Rate</th>
                            <th>Avg Count/Gene</th>
                        </tr>
                    </thead>
                    <tbody>
"""
    
    if count_analysis:
        total_genes = count_analysis['gene_stats']['total_genes']
        for sample in samples:
            data = count_analysis['samples'].get(sample, {})
            total_counts = data.get('total_counts', 0)
            genes_detected = data.get('genes_detected', 0)
            detection_rate = (genes_detected / total_genes * 100) if total_genes > 0 else 0
            avg_count = (total_counts / total_genes) if total_genes > 0 else 0
            
            html += f"""
                        <tr>
                            <td><strong>{sample}</strong></td>
                            <td>{total_counts:,}</td>
                            <td>{genes_detected:,}</td>
                            <td>
                                {detection_rate:.1f}%
                                <div class="progress-bar">
                                    <div class="progress-fill" style="width: {min(detection_rate, 100)}%">{detection_rate:.1f}%</div>
                                </div>
                            </td>
                            <td>{avg_count:.1f}</td>
                        </tr>
"""
    
    html += """
                    </tbody>
                </table>
            </div>
            
            <!-- 고발현 유전자 -->
            <div class="section">
                <h2 class="section-title">🔝 Top Expressed Genes</h2>
"""
    
    if count_analysis:
        for sample in samples:
            data = count_analysis['samples'].get(sample, {})
            high_expr = data.get('high_expression', [])
            
            if high_expr:
                html += f"""
                <h3 style="color: #667eea; margin-top: 20px;">{sample}</h3>
                <div class="gene-list">
"""
                for gene_id, count in high_expr[:10]:
                    html += f"""
                    <div class="gene-item">
                        <span><strong>{gene_id}</strong></span>
                        <span class="metric-good">{count:,} reads</span>
                    </div>
"""
                html += """
                </div>
"""
    
    html += """
            </div>
            
            <!-- 파일 크기 -->
            <div class="section">
                <h2 class="section-title">💾 File Sizes</h2>
                <table>
                    <thead>
                        <tr>
                            <th>Sample</th>
                            <th>Raw R1</th>
                            <th>Raw R2</th>
                            <th>Trimmed R1</th>
                            <th>Trimmed R2</th>
                            <th>Aligned BAM</th>
                        </tr>
                    </thead>
                    <tbody>
"""
    
    for sample in samples:
        raw_r1 = RAW_FILES.get(sample, {}).get('r1', '')
        raw_r2 = RAW_FILES.get(sample, {}).get('r2', '')
        trimmed_r1 = os.path.join(TRIMMED_DIR, f'{sample}_1.fastq.gz')
        trimmed_r2 = os.path.join(TRIMMED_DIR, f'{sample}_2.fastq.gz')
        bam = os.path.join(ALIGNED_DIR, sample, 'Aligned.sortedByCoord.out.bam')

        html += f"""
                        <tr>
                            <td><strong>{sample}</strong></td>
                            <td>{get_file_size(raw_r1)}</td>
                            <td>{get_file_size(raw_r2)}</td>
                            <td>{get_file_size(trimmed_r1)}</td>
                            <td>{get_file_size(trimmed_r2)}</td>
                            <td>{get_file_size(bam)}</td>
                        </tr>
"""

    counts_matrix = os.path.join(COUNTS_DIR, 'counts_matrix.txt')
    html += f"""
                    </tbody>
                </table>
                <div style="margin-top: 20px;">
                    <strong>Count Matrix:</strong> {get_file_size(counts_matrix)}
                </div>
            </div>
"""

    # ── Pipeline QC Evaluation Summary (English) ──────────────────
    html += _render_pipeline_qc_section(PIPELINE_QC_JSON)

    html += f"""
        </div>

        <div class="footer">
            <p>Generated by RNA-seq Pipeline QC Report Generator</p>
            <p>Pipeline Version: 1.0 | Report Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        </div>
    </div>
</body>
</html>
"""
    
    return html

def main():
    """메인 실행 함수"""
    # 로깅 설정
    if LOG_FILE and SCRIPT_MODE:
        import sys
        log_f = open(LOG_FILE, 'w')
        sys.stdout = log_f
        sys.stderr = log_f
    
    print("🔬 RNA-seq Pipeline QC Report Generator")
    print("=" * 60)
    print(f"Script mode: {SCRIPT_MODE}")
    print(f"RESULTS_DIR: {RESULTS_DIR}")
    print(f"LOGS_DIR: {LOGS_DIR}")
    print(f"BASE_DIR: {BASE_DIR}")
    
    # 샘플 이름 수집
    samples = get_sample_names()
    
    if not samples:
        print(f"❌ No samples found in {os.path.join(RESULTS_DIR, 'trimmed')}")
        print("   Make sure the pipeline has been run successfully.")
        if LOG_FILE and SCRIPT_MODE:
            log_f.close()
        return
    
    print(f"✓ Found {len(samples)} sample(s): {', '.join(samples)}")
    print(f"✓ Top genes to show: {TOP_GENES}")
    
    # HTML 리포트 생성
    print("\n📝 Generating HTML report...")
    html_content = generate_html_report(samples)
    
    # 파일 저장
    with open(OUTPUT_FILE, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    print(f"✅ Report saved: {OUTPUT_FILE}")
    
    if not SCRIPT_MODE:
        print(f"\n🌐 Open the report in your browser:")
        print(f"   file://{os.path.abspath(OUTPUT_FILE)}")
    
    # 터미널 요약도 출력
    print("\n" + "=" * 60)
    print("📊 Quick Summary")
    print("=" * 60)
    
    for sample in samples:
        print(f"\n{sample}:")
        star_data = parse_star_log(sample)
        if star_data:
            print(f"  Mapping rate: {star_data.get('unique_mapped_pct', 'N/A')}")
        
        count_analysis = analyze_count_matrix()
        if count_analysis and sample in count_analysis['samples']:
            data = count_analysis['samples'][sample]
            print(f"  Genes detected: {data['genes_detected']:,}")
            print(f"  Total counts: {data['total_counts']:,}")
    
    if LOG_FILE and SCRIPT_MODE:
        log_f.close()

if __name__ == '__main__':
    main()
