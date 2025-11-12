import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from scipy import stats
import io
import json
from datetime import datetime
from typing import Dict, List, Optional

# ==================== CONFIG ====================
st.set_page_config(page_title="qPCR Analysis Pro", layout="wide", initial_sidebar_state="expanded")

# Session state initialization
for key in ['data', 'processed_data', 'sample_mapping', 'sample_merge_map', 'analysis_templates', 
            'graphs', 'excluded_wells', 'excluded_samples', 'selected_efficacy', 'hk_gene', 
            'included_samples_for_analysis', 'graph_settings']:
    if key not in st.session_state:
        st.session_state[key] = {} if key in ['sample_mapping', 'sample_merge_map', 'analysis_templates', 
                                               'graphs', 'graph_settings', 'included_samples_for_analysis', 
                                               'processed_data'] else (set() if 'excluded' in key else None)

# Efficacy configurations
EFFICACY_CONFIG = {
    '탄력': {'genes': ['COL1A1', 'ELN', 'FBN-1', 'FBN1'], 'cell': 'HS68', 
            'controls': {'baseline': 'Non-treated', 'positive': 'TGFb'}},
    '항노화': {'genes': ['COL1A1', 'COL1', 'MMP-1', 'MMP1'], 'cell': 'HS68',
              'controls': {'baseline': 'Non-treated', 'negative': 'UVB', 'positive': 'UVB+TGFb'},
              'expected': {'COL1A1': 'up', 'MMP-1': 'down'}},
    '보습': {'genes': ['AQP3', 'HAS3'], 'cell': 'HaCaT',
           'controls': {'baseline': 'Non-treated', 'positive': 'Retinoic acid'}},
    '장벽': {'genes': ['FLG', 'CLDN', 'IVL'], 'cell': 'HaCaT',
           'controls': {'baseline': 'Non-treated', 'positive': 'Retinoic acid'}},
    '표피증식': {'genes': ['KI67', 'PCNA'], 'cell': 'HaCaT',
              'controls': {'baseline': 'Non-treated', 'positive': 'TGFb/FBS'}},
    '멜라닌억제': {'genes': ['MITF', 'TYR'], 'cell': 'B16F10',
                'controls': {'baseline': 'Non-treated', 'negative': 'α-MSH', 'positive': 'α-MSH+Arbutin'},
                'expected': {'MITF': 'down', 'TYR': 'down'}},
    '진정': {'genes': ['IL1B', 'IL-1β', 'IL6', 'TNFA', 'TNFα'], 'cell': 'HaCaT',
           'controls': {'baseline': 'Non-treated', 'negative': 'Inflammation', 'positive': 'Inflam+Dex'},
           'expected': {'IL1B': 'down', 'IL6': 'down', 'TNFA': 'down'}},
    '지질억제': {'genes': ['SREBPA', 'SREBPa', 'SREBPC', 'SREBPc'], 'cell': 'SZ95',
              'controls': {'baseline': 'Non-treated', 'negative': 'IGF', 'positive': 'IGF+inhibitor'},
              'expected': {'SREBPA': 'down', 'SREBPC': 'down'}},
    '냉감': {'genes': ['TRPM8', 'CIRBP'], 'cell': 'HaCaT',
           'controls': {'baseline': 'Non-treated', 'positive': 'Menthol'}}
}

# ==================== PARSER ====================
class QPCRParser:
    @staticmethod
    def detect_format(df):
        for idx, row in df.iterrows():
            row_str = ' '.join(row.astype(str).values)
            if 'Well Position' in row_str:
                return 'format1', idx
            elif row.iloc[0] == 'Well' and 'Sample Name' in row_str:
                return ('format2' if 'Cт' in row_str else 'format1'), idx
        return 'unknown', 0
    
    @staticmethod
    def parse_format(df, start, ct_col_name='CT'):
        df = df.iloc[start:].reset_index(drop=True)
        df.columns = df.iloc[0]
        df = df.iloc[1:].reset_index(drop=True)
        
        well_col = next((c for c in ['Well Position', 'Well'] if c in df.columns), df.columns[0])
        ct_col = next((c for c in ['CT', 'Ct', 'Cт'] if c in df.columns), None)
        if not ct_col:
            return None
        
        return pd.DataFrame({
            'Well': df[well_col],
            'Sample': df.get('Sample Name', df.iloc[:, 2]),
            'Target': df.get('Target Name', df.iloc[:, 3]),
            'CT': pd.to_numeric(df[ct_col], errors='coerce')
        }).dropna(subset=['CT']).query('Sample.notna() & Target.notna()')
    
    @staticmethod
    def parse(file):
        try:
            df = None
            for enc in ['utf-8', 'latin-1', 'cp1252']:
                try:
                    df = pd.read_csv(file, encoding=enc, low_memory=False, skip_blank_lines=False)
                    break
                except:
                    continue
            if df is None:
                return None
            fmt, start = QPCRParser.detect_format(df)
            return QPCRParser.parse_format(df, start) if fmt != 'unknown' else None
        except Exception as e:
            st.error(f"Parse error: {e}")
            return None

# ==================== ANALYSIS ENGINE ====================
class AnalysisEngine:
    @staticmethod
    def calculate_ddct(data: pd.DataFrame, hk_gene: str, baseline_condition: str,
                       sample_mapping: dict, included_samples: Dict[str, List[str]]) -> Dict[str, pd.DataFrame]:
        
        data = data.copy()
        data['Condition'] = data['Sample'].map(lambda x: sample_mapping.get(x, {}).get('condition', x))
        data['Group'] = data['Sample'].map(lambda x: sample_mapping.get(x, {}).get('group', 'Treatment'))
        
        gene_results = {}
        
        for target in data['Target'].unique():
            if target.upper() in [hk_gene.upper(), 'ACTIN', 'B-ACTIN', 'GAPDH', 'ACTB']:
                continue
            
            included = included_samples.get(target, [])
            if not included or baseline_condition not in included:
                continue
            
            target_data = data[(data['Target'] == target) & (data['Condition'].isin(included))]
            
            # Get baseline ΔCt
            baseline_target = target_data[target_data['Condition'] == baseline_condition]
            baseline_hk = data[(data['Condition'] == baseline_condition) & (data['Target'] == hk_gene)]
            if len(baseline_target) == 0 or len(baseline_hk) == 0:
                continue
            
            baseline_delta_ct = baseline_target['CT'].mean() - baseline_hk['CT'].mean()
            
            results = []
            for condition in target_data['Condition'].unique():
                cond_target = target_data[target_data['Condition'] == condition]
                cond_hk = data[(data['Condition'] == condition) & (data['Target'] == hk_gene)]
                
                if len(cond_hk) == 0:
                    continue
                
                target_ct_values = cond_target['CT'].values
                delta_ct = target_ct_values.mean() - cond_hk['CT'].mean()
                ddct = delta_ct - baseline_delta_ct
                rel_expr = 2 ** (-ddct)
                
                sem = target_ct_values.std() / np.sqrt(len(target_ct_values)) if len(target_ct_values) > 1 else 0
                
                results.append({
                    'Target': target, 'Condition': condition,
                    'Group': cond_target['Group'].iloc[0],
                    'n_replicates': len(target_ct_values),
                    'Target_Ct_Values': target_ct_values,
                    'Delta_Ct': delta_ct, 'Delta_Delta_Ct': ddct,
                    'Fold_Change': rel_expr, 'SEM': sem
                })
            
            if results:
                gene_results[target] = pd.DataFrame(results)
        
        return gene_results
    
    @staticmethod
    def calculate_statistics(gene_data: pd.DataFrame, ref_condition: str) -> pd.DataFrame:
        results = gene_data.copy()
        results['p_value'] = np.nan
        results['significance'] = ''
        
        ref_data = results[results['Condition'] == ref_condition]
        if len(ref_data) == 0:
            return results
        
        ref_values = ref_data.iloc[0]['Target_Ct_Values']
        
        for idx, row in results.iterrows():
            if row['Condition'] == ref_condition:
                results.at[idx, 'p_value'] = 1.0
                continue
            
            sample_values = row['Target_Ct_Values']
            
            try:
                # Two-tailed t-test on raw CT values
                if len(ref_values) > 1 and len(sample_values) > 1:
                    _, p_val = stats.ttest_ind(ref_values, sample_values)
                elif len(ref_values) == 1 and len(sample_values) > 1:
                    _, p_val = stats.ttest_1samp(sample_values, ref_values[0])
                elif len(sample_values) == 1 and len(ref_values) > 1:
                    _, p_val = stats.ttest_1samp(ref_values, sample_values[0])
                else:
                    p_val = np.nan
                
                results.at[idx, 'p_value'] = p_val
                
                if p_val < 0.001:
                    results.at[idx, 'significance'] = '***'
                elif p_val < 0.01:
                    results.at[idx, 'significance'] = '**'
                elif p_val < 0.05:
                    results.at[idx, 'significance'] = '*'
            except:
                results.at[idx, 'p_value'] = np.nan
        
        return results

# ==================== GRAPH GENERATOR ====================
class GraphGenerator:
    @staticmethod
    def create_graph(data: pd.DataFrame, gene: str, settings: dict) -> go.Figure:
        gene_data = data[data['Condition'].isin(settings['included_conditions'])].copy()
        
        if settings.get('custom_order'):
            order_map = {c: i for i, c in enumerate(settings['custom_order'])}
            gene_data['sort_key'] = gene_data['Condition'].map(lambda x: order_map.get(x, 999))
            gene_data = gene_data.sort_values('sort_key')
        
        colors = [settings['condition_colors'].get(c, '#636EFA') for c in gene_data['Condition']]
        
        # Create significance text with custom formatting
        sig_text = []
        for _, row in gene_data.iterrows():
            if settings.get('show_significance', True) and row['significance']:
                sig_text.append(row['significance'])
            else:
                sig_text.append('')
        
        fig = go.Figure()
        
        fig.add_trace(go.Bar(
            x=gene_data['Condition'],
            y=gene_data['Fold_Change'],
            error_y=dict(
                type='data',
                array=gene_data['SEM'] * 1.96,
                visible=settings.get('show_error', True),
                thickness=settings.get('error_thickness', 2),
                width=settings.get('error_width', 4),
                color=settings.get('error_color', 'black')
            ),
            text=sig_text,
            textposition=settings.get('sig_position', 'outside'),
            textfont=dict(
                size=settings.get('sig_size', 16),
                color=settings.get('sig_color', 'black'),
                family=settings.get('sig_font', 'Arial')
            ),
            marker=dict(
                color=colors,
                line=dict(
                    color=settings.get('bar_border_color', 'black'),
                    width=settings.get('bar_border_width', 0)
                ),
                opacity=settings.get('bar_opacity', 1.0)
            ),
            width=settings.get('bar_width', 0.8),
            showlegend=False
        ))
        
        # Reference line
        if settings.get('show_reference_line', True):
            fig.add_hline(
                y=settings.get('reference_value', 1.0),
                line_dash="dash",
                line_color=settings.get('reference_line_color', 'gray'),
                line_width=settings.get('reference_line_width', 2)
            )
        
        # Layout
        fig.update_layout(
            title=dict(
                text=settings.get('title', f"{gene} Expression"),
                font=dict(
                    size=settings.get('title_size', 20),
                    color=settings.get('title_color', 'black'),
                    family=settings.get('title_font', 'Arial')
                ),
                x=settings.get('title_x', 0.5),
                xanchor='center'
            ),
            xaxis=dict(
                title=dict(
                    text=settings.get('xlabel', 'Condition'),
                    font=dict(size=settings.get('axis_label_size', 14))
                ),
                showgrid=settings.get('show_x_grid', False),
                gridcolor=settings.get('grid_color', 'lightgray'),
                tickangle=settings.get('x_tick_angle', -45),
                tickfont=dict(size=settings.get('tick_size', 12))
            ),
            yaxis=dict(
                title=dict(
                    text=settings.get('ylabel', 'Fold Change'),
                    font=dict(size=settings.get('axis_label_size', 14))
                ),
                showgrid=settings.get('show_y_grid', True),
                gridcolor=settings.get('grid_color', 'lightgray'),
                tickfont=dict(size=settings.get('tick_size', 12)),
                range=settings.get('y_range')
            ),
            template=settings.get('template', 'plotly_white'),
            height=settings.get('height', 600),
            width=settings.get('width', 1000),
            plot_bgcolor=settings.get('plot_bgcolor', 'white'),
            paper_bgcolor=settings.get('paper_bgcolor', 'white'),
            margin=dict(l=80, r=80, t=100, b=100)
        )
        
        return fig

# ==================== EXPORT ====================
def export_excel(raw_data: pd.DataFrame, processed_data: Dict, params: dict, mapping: dict) -> bytes:
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
        pd.DataFrame([params]).to_excel(writer, sheet_name='Parameters', index=False)
        pd.DataFrame([{'Original': k, **v} for k, v in mapping.items()]).to_excel(
            writer, sheet_name='Mapping', index=False)
        raw_data.to_excel(writer, sheet_name='Raw_Data', index=False)
        
        for gene, df in processed_data.items():
            export_df = df.drop(columns=['Target_Ct_Values'], errors='ignore')
            export_df.to_excel(writer, sheet_name=f"{gene}"[:31], index=False)
    
    return output.getvalue()

# ==================== UI ====================
st.title("🧬 qPCR Analysis Pro")

# Sidebar
with st.sidebar:
    st.header("📖 Guide")
    st.markdown("1. Upload & merge samples\n2. Map conditions\n3. Analyze\n4. Customize graphs\n5. Export")
    
    if st.button("💾 Save Template") and st.session_state.sample_mapping:
        name = st.text_input("Name:", key="save_name")
        if name:
            st.session_state.analysis_templates[name] = {
                'mapping': st.session_state.sample_mapping.copy(),
                'settings': st.session_state.graph_settings.copy()
            }
            st.success(f"Saved '{name}'")

# Tabs
tab1, tab2, tab3, tab4, tab5 = st.tabs(["📁 Data", "🗺️ Mapping", "🔬 Analysis", "📊 Graphs", "📤 Export"])

# ==================== TAB 1: DATA ====================
with tab1:
    st.header("Upload & Filter Data")
    
    files = st.file_uploader("Upload CSV files", type=['csv'], accept_multiple_files=True)
    
    if files:
        all_data = []
        for file in files:
            parsed = QPCRParser.parse(file)
            if parsed is not None:
                parsed['Source'] = file.name
                all_data.append(parsed)
                st.success(f"✅ {file.name}: {len(parsed)} wells")
        
        if all_data:
            st.session_state.data = pd.concat(all_data, ignore_index=True)
            
            # Sample merger
            st.subheader("🔗 Merge Samples")
            st.markdown("*Combine differently-labeled samples from multiple files*")
            
            all_samples = sorted(st.session_state.data['Sample'].unique())
            
            if st.checkbox("Enable sample merging"):
                merge_groups = st.number_input("Number of merge groups", 1, 10, 1)
                
                for i in range(merge_groups):
                    with st.expander(f"Merge Group {i+1}"):
                        samples_to_merge = st.multiselect(
                            f"Select samples to merge into one:",
                            all_samples,
                            key=f"merge_{i}"
                        )
                        if samples_to_merge:
                            new_name = st.text_input(f"Merged name:", samples_to_merge[0], key=f"mergename_{i}")
                            if st.button(f"Apply merge {i+1}", key=f"apply_{i}"):
                                for sample in samples_to_merge:
                                    st.session_state.sample_merge_map[sample] = new_name
                                st.success(f"Merged {len(samples_to_merge)} samples → '{new_name}'")
                
                # Apply merging
                if st.session_state.sample_merge_map:
                    st.session_state.data['Sample'] = st.session_state.data['Sample'].map(
                        lambda x: st.session_state.sample_merge_map.get(x, x)
                    )
            
            # Housekeeping gene
            hk_genes = [g for g in st.session_state.data['Target'].unique() 
                       if g.upper() in ['ACTIN', 'B-ACTIN', 'GAPDH', 'ACTB']]
            if hk_genes:
                st.session_state.hk_gene = st.selectbox("🔬 Housekeeping Gene", hk_genes)
            
            # Filters
            st.subheader("🎯 Filters")
            
            col1, col2 = st.columns(2)
            with col1:
                samples = sorted(st.session_state.data['Sample'].unique())
                selected_samples = st.multiselect("Include samples:", samples, default=samples)
                st.session_state.excluded_samples = set(samples) - set(selected_samples)
            
            with col2:
                genes = [g for g in sorted(st.session_state.data['Target'].unique()) 
                        if g.upper() not in ['ACTIN', 'B-ACTIN', 'GAPDH', 'ACTB']]
                selected_genes = st.multiselect("Include genes:", genes, default=genes)
            
            # Show all filtered data
            st.subheader("📊 All Filtered Data")
            
            display_data = st.session_state.data[
                st.session_state.data['Sample'].isin(selected_samples) &
                st.session_state.data['Target'].isin(selected_genes + [st.session_state.hk_gene])
            ]
            
            # Group by sample and show statistics
            st.markdown("**Data by Sample:**")
            for sample in selected_samples:
                with st.expander(f"📍 {sample}"):
                    sample_data = display_data[display_data['Sample'] == sample]
                    
                    # Show CT values for each gene
                    summary = sample_data.groupby('Target')['CT'].agg(['mean', 'std', 'count']).round(2)
                    st.dataframe(summary, width='stretch')
                    
                    # Show raw data
                    st.caption("Raw data:")
                    st.dataframe(sample_data[['Well', 'Target', 'CT', 'Source']], width='stretch')

# ==================== TAB 2: MAPPING ====================
with tab2:
    st.header("Sample Mapping")
    
    if st.session_state.data is not None:
        # Efficacy selection
        detected = set(st.session_state.data['Target'].unique())
        suggested = next((e for e, c in EFFICACY_CONFIG.items() 
                         if any(g in detected for g in c['genes'])), None)
        
        efficacy = st.selectbox("🎯 Efficacy Type", list(EFFICACY_CONFIG.keys()),
                               index=list(EFFICACY_CONFIG.keys()).index(suggested) if suggested else 0)
        st.session_state.selected_efficacy = efficacy
        
        # Mapping table with headers
        st.subheader("🗺️ Condition Mapping")
        
        samples = [s for s in sorted(st.session_state.data['Sample'].unique()) 
                  if s not in st.session_state.excluded_samples]
        
        # Initialize mappings
        for sample in samples:
            if sample not in st.session_state.sample_mapping:
                st.session_state.sample_mapping[sample] = {
                    'condition': sample, 'group': 'Treatment', 'concentration': ''
                }
        
        # Display as table with headers
        st.markdown("**Edit sample information:**")
        
        # Headers
        col0, col1, col2, col3 = st.columns([1, 3, 2, 2])
        col0.markdown("**Original Sample**")
        col1.markdown("**Condition Name**")
        col2.markdown("**Group Type**")
        col3.markdown("**Concentration**")
        
        # Data rows
        group_types = ['Baseline', 'Negative Control', 'Positive Control', 'Treatment']
        
        for sample in samples:
            col0, col1, col2, col3 = st.columns([1, 3, 2, 2])
            
            col0.text(sample)
            
            cond = col1.text_input(
                "Condition",
                st.session_state.sample_mapping[sample]['condition'],
                key=f"cond_{sample}",
                label_visibility="collapsed"
            )
            st.session_state.sample_mapping[sample]['condition'] = cond
            
            grp = col2.selectbox(
                "Group",
                group_types,
                index=group_types.index(st.session_state.sample_mapping[sample]['group']) 
                      if st.session_state.sample_mapping[sample]['group'] in group_types else 0,
                key=f"grp_{sample}",
                label_visibility="collapsed"
            )
            st.session_state.sample_mapping[sample]['group'] = grp
            
            conc = col3.text_input(
                "Conc",
                st.session_state.sample_mapping[sample]['concentration'],
                key=f"conc_{sample}",
                label_visibility="collapsed"
            )
            st.session_state.sample_mapping[sample]['concentration'] = conc

# ==================== TAB 3: ANALYSIS ====================
with tab3:
    st.header("Analysis")
    
    if st.session_state.data is not None and st.session_state.sample_mapping and st.session_state.hk_gene:
        
        col1, col2 = st.columns(2)
        
        with col1:
            baseline_samples = [k for k, v in st.session_state.sample_mapping.items() 
                              if v['group'] == 'Baseline']
            baseline_sample = st.selectbox("🎯 Baseline (Fold Change = 1.0)", 
                                          baseline_samples if baseline_samples else 
                                          list(st.session_state.sample_mapping.keys()))
            baseline_condition = st.session_state.sample_mapping[baseline_sample]['condition']
        
        with col2:
            all_conditions = [v['condition'] for v in st.session_state.sample_mapping.values()]
            pvalue_condition = st.selectbox("📊 P-value Reference", all_conditions,
                                           index=all_conditions.index(baseline_condition) 
                                           if baseline_condition in all_conditions else 0)
        
        # Per-gene sample selection
        st.subheader("🧬 Select Samples Per Gene")
        
        genes = [g for g in st.session_state.data['Target'].unique() 
                if g.upper() not in ['ACTIN', 'B-ACTIN', 'GAPDH', 'ACTB']]
        
        for gene in genes:
            if gene not in st.session_state.included_samples_for_analysis:
                st.session_state.included_samples_for_analysis[gene] = all_conditions.copy()
            
            selected = st.multiselect(
                f"📍 {gene}:",
                all_conditions,
                default=st.session_state.included_samples_for_analysis[gene],
                key=f"incl_{gene}"
            )
            st.session_state.included_samples_for_analysis[gene] = selected
        
        # Run analysis
        if st.button("🔬 Run Analysis", type="primary"):
            with st.spinner("Analyzing..."):
                gene_results = AnalysisEngine.calculate_ddct(
                    st.session_state.data,
                    st.session_state.hk_gene,
                    baseline_condition,
                    st.session_state.sample_mapping,
                    st.session_state.included_samples_for_analysis
                )
                
                if gene_results:
                    for gene, df in gene_results.items():
                        gene_results[gene] = AnalysisEngine.calculate_statistics(df, pvalue_condition)
                    
                    st.session_state.processed_data = gene_results
                    st.success(f"✅ Analyzed {len(gene_results)} genes")
        
        # Display results
        if st.session_state.processed_data:
            st.subheader("📊 Results")
            
            for gene, df in st.session_state.processed_data.items():
                with st.expander(f"🧬 {gene}", expanded=False):
                    display_df = df.drop(columns=['Target_Ct_Values'], errors='ignore')
                    display_cols = ['Condition', 'Group', 'Fold_Change', 'p_value', 
                                  'significance', 'n_replicates', 'SEM']
                    
                    styled = display_df[[c for c in display_cols if c in display_df.columns]].style\
                        .background_gradient(subset=['Fold_Change'], cmap='RdYlGn', vmin=0, vmax=3)\
                        .format({'Fold_Change': '{:.3f}', 'p_value': '{:.4f}', 'SEM': '{:.3f}'})
                    
                    st.dataframe(styled, width='stretch')

# ==================== TAB 4: GRAPHS ====================
with tab4:
    st.header("Interactive Graph Editor")
    
    if st.session_state.processed_data:
        
        # Initialize settings
        for gene in st.session_state.processed_data.keys():
            if gene not in st.session_state.graph_settings:
                colors = px.colors.qualitative.Plotly
                idx = list(st.session_state.processed_data.keys()).index(gene)
                
                st.session_state.graph_settings[gene] = {
                    'title': f"{gene} Expression", 'xlabel': 'Condition',
                    'ylabel': 'Fold Change', 'title_size': 20, 'axis_label_size': 14,
                    'tick_size': 12, 'sig_size': 16, 'width': 1000, 'height': 600,
                    'template': 'plotly_white', 'show_error': True, 'show_significance': True,
                    'show_x_grid': False, 'show_y_grid': True, 'bar_width': 0.8,
                    'show_reference_line': True, 'reference_value': 1.0,
                    'condition_colors': {}, 'included_conditions': [],
                    'custom_order': None, 'x_tick_angle': -45, 'sig_position': 'outside',
                    'sig_color': 'black', 'sig_font': 'Arial', 'bar_border_width': 0,
                    'bar_border_color': 'black', 'bar_opacity': 1.0, 'error_thickness': 2,
                    'error_width': 4, 'error_color': 'black', 'grid_color': 'lightgray',
                    'plot_bgcolor': 'white', 'paper_bgcolor': 'white', 'y_range': None,
                    'title_color': 'black', 'title_font': 'Arial', 'title_x': 0.5,
                    'reference_line_color': 'gray', 'reference_line_width': 2
                }
                
                # Initialize included conditions
                st.session_state.graph_settings[gene]['included_conditions'] = \
                    list(st.session_state.processed_data[gene]['Condition'].unique())
        
        # Quick presets
        st.subheader("⚡ Quick Presets")
        col1, col2, col3, col4 = st.columns(4)
        
        if col1.button("📄 Publication"):
            for gene in st.session_state.graph_settings:
                st.session_state.graph_settings[gene].update({
                    'template': 'simple_white', 'width': 1000, 'height': 600,
                    'bar_border_width': 1, 'title_size': 18, 'axis_label_size': 14
                })
            st.rerun()
        
        if col2.button("🎨 Presentation"):
            for gene in st.session_state.graph_settings:
                st.session_state.graph_settings[gene].update({
                    'template': 'presentation', 'width': 1200, 'height': 700,
                    'title_size': 24, 'axis_label_size': 18, 'tick_size': 14
                })
            st.rerun()
        
        if col3.button("🌙 Dark"):
            for gene in st.session_state.graph_settings:
                st.session_state.graph_settings[gene].update({
                    'template': 'plotly_dark', 'plot_bgcolor': '#111', 'paper_bgcolor': '#111',
                    'title_color': 'white', 'sig_color': 'white'
                })
            st.rerun()
        
        if col4.button("🔄 Reset"):
            st.session_state.graph_settings = {}
            st.rerun()
        
        st.markdown("---")
        
        # Gene-specific customization
        for gene in st.session_state.processed_data.keys():
            st.markdown(f"## 🧬 {gene}")
            
            settings = st.session_state.graph_settings[gene]
            gene_data = st.session_state.processed_data[gene]
            
            # Tabs for organization
            tab_samples, tab_bars, tab_text, tab_style, tab_colors = st.tabs([
                "📋 Samples", "📊 Bars", "📝 Text", "🎨 Style", "🌈 Colors"
            ])
            
            # TAB: Samples
            with tab_samples:
                all_conds = list(gene_data['Condition'].unique())
                settings['included_conditions'] = st.multiselect(
                    "Display these conditions:",
                    all_conds,
                    default=settings['included_conditions'],
                    key=f"incl_{gene}"
                )
                
                if st.checkbox(f"Custom order", key=f"ord_{gene}"):
                    settings['custom_order'] = st.multiselect(
                        "Drag to reorder:",
                        settings['included_conditions'],
                        default=settings['included_conditions'],
                        key=f"order_{gene}"
                    )
            
            # TAB: Bars
            with tab_bars:
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    st.markdown("**Bar Appearance**")
                    settings['bar_width'] = st.slider("Width", 0.1, 1.5, 
                                                     float(settings.get('bar_width', 0.8)),
                                                     key=f"barw_{gene}")
                    settings['bar_opacity'] = st.slider("Opacity", 0.0, 1.0,
                                                       float(settings.get('bar_opacity', 1.0)),
                                                       key=f"barop_{gene}")
                    settings['bar_border_width'] = st.slider("Border Width", 0, 5,
                                                            int(settings.get('bar_border_width', 0)),
                                                            key=f"barbord_{gene}")
                
                with col2:
                    st.markdown("**Error Bars**")
                    settings['show_error'] = st.checkbox("Show Error Bars",
                                                        bool(settings.get('show_error', True)),
                                                        key=f"err_{gene}")
                    if settings['show_error']:
                        settings['error_thickness'] = st.slider("Thickness", 1, 5,
                                                               int(settings.get('error_thickness', 2)),
                                                               key=f"errthick_{gene}")
                        settings['error_width'] = st.slider("Cap Width", 2, 10,
                                                           int(settings.get('error_width', 4)),
                                                           key=f"errw_{gene}")
                
                with col3:
                    st.markdown("**Reference Line**")
                    settings['show_reference_line'] = st.checkbox("Show Reference",
                                                                 bool(settings.get('show_reference_line', True)),
                                                                 key=f"ref_{gene}")
                    if settings['show_reference_line']:
                        settings['reference_value'] = st.number_input("Value",
                                                                      value=float(settings.get('reference_value', 1.0)),
                                                                      key=f"refv_{gene}")
                        settings['reference_line_width'] = st.slider("Line Width", 1, 5,
                                                                     int(settings.get('reference_line_width', 2)),
                                                                     key=f"refw_{gene}")
            
            # TAB: Text
            with tab_text:
                col1, col2 = st.columns(2)
                
                with col1:
                    st.markdown("**Labels**")
                    settings['title'] = st.text_input("Title",
                                                     str(settings.get('title', f"{gene} Expression")),
                                                     key=f"title_{gene}")
                    settings['xlabel'] = st.text_input("X-axis",
                                                      str(settings.get('xlabel', 'Condition')),
                                                      key=f"xlabel_{gene}")
                    settings['ylabel'] = st.text_input("Y-axis",
                                                      str(settings.get('ylabel', 'Fold Change')),
                                                      key=f"ylabel_{gene}")
                
                with col2:
                    st.markdown("**Sizes**")
                    settings['title_size'] = st.slider("Title Size", 12, 36,
                                                      int(settings.get('title_size', 20)),
                                                      key=f"titsize_{gene}")
                    settings['axis_label_size'] = st.slider("Axis Label Size", 10, 24,
                                                           int(settings.get('axis_label_size', 14)),
                                                           key=f"axsize_{gene}")
                    settings['tick_size'] = st.slider("Tick Size", 8, 20,
                                                     int(settings.get('tick_size', 12)),
                                                     key=f"ticksize_{gene}")
                    settings['x_tick_angle'] = st.slider("X Tick Angle", -90, 90,
                                                         int(settings.get('x_tick_angle', -45)),
                                                         key=f"tickang_{gene}")
                
                st.markdown("**Significance Stars**")
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    settings['show_significance'] = st.checkbox("Show Stars",
                                                               bool(settings.get('show_significance', True)),
                                                               key=f"sig_{gene}")
                
                if settings['show_significance']:
                    with col2:
                        positions = ['outside', 'inside', 'auto']
                        curr_pos = settings.get('sig_position', 'outside')
                        settings['sig_position'] = st.selectbox("Position", positions,
                                                               index=positions.index(curr_pos) if curr_pos in positions else 0,
                                                               key=f"sigpos_{gene}")
                    
                    with col3:
                        settings['sig_size'] = st.slider("Star Size", 10, 30,
                                                        int(settings.get('sig_size', 16)),
                                                        key=f"sigsize_{gene}")
            
            # TAB: Style
            with tab_style:
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    st.markdown("**Theme**")
                    templates = ['plotly_white', 'plotly', 'plotly_dark', 'seaborn', 
                                'simple_white', 'presentation', 'ggplot2', 'none']
                    curr_template = settings.get('template', 'plotly_white')
                    settings['template'] = st.selectbox("Template", templates,
                                                       index=templates.index(curr_template) if curr_template in templates else 0,
                                                       key=f"temp_{gene}")
                    
                    settings['width'] = st.number_input("Width (px)", 400, 2000,
                                                       int(settings.get('width', 1000)),
                                                       key=f"w_{gene}")
                    settings['height'] = st.number_input("Height (px)", 300, 1500,
                                                        int(settings.get('height', 600)),
                                                        key=f"h_{gene}")
                
                with col2:
                    st.markdown("**Grid**")
                    settings['show_x_grid'] = st.checkbox("X Grid",
                                                         bool(settings.get('show_x_grid', False)),
                                                         key=f"xgrid_{gene}")
                    settings['show_y_grid'] = st.checkbox("Y Grid",
                                                         bool(settings.get('show_y_grid', True)),
                                                         key=f"ygrid_{gene}")
                    settings['grid_color'] = st.color_picker("Grid Color",
                                                            str(settings.get('grid_color', 'lightgray')),
                                                            key=f"gridc_{gene}")
                
                with col3:
                    st.markdown("**Y-axis Range**")
                    use_custom = st.checkbox("Custom Range",
                                           value=settings.get('y_range') is not None,
                                           key=f"yrange_{gene}")
                    if use_custom:
                        curr_range = settings.get('y_range', [0, 3])
                        y_min = st.number_input("Y Min", value=float(curr_range[0]), key=f"ymin_{gene}")
                        y_max = st.number_input("Y Max", value=float(curr_range[1]), key=f"ymax_{gene}")
                        settings['y_range'] = [y_min, y_max]
                    else:
                        settings['y_range'] = None
            
            # TAB: Colors
            with tab_colors:
                st.markdown("**Individual Bar Colors**")
                
                # Initialize colors
                if not settings['condition_colors']:
                    default_colors = px.colors.qualitative.Plotly
                    for i, cond in enumerate(settings['included_conditions']):
                        settings['condition_colors'][cond] = default_colors[i % len(default_colors)]
                
                cols = st.columns(3)
                for i, cond in enumerate(settings['included_conditions']):
                    with cols[i % 3]:
                        if cond not in settings['condition_colors']:
                            settings['condition_colors'][cond] = '#636EFA'
                        
                        settings['condition_colors'][cond] = st.color_picker(
                            cond,
                            settings['condition_colors'][cond],
                            key=f"col_{gene}_{cond}"
                        )
                
                st.markdown("**Other Colors**")
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    settings['title_color'] = st.color_picker("Title",
                                                             str(settings.get('title_color', 'black')),
                                                             key=f"titc_{gene}")
                    settings['sig_color'] = st.color_picker("Significance",
                                                           str(settings.get('sig_color', 'black')),
                                                           key=f"sigc_{gene}")
                
                with col2:
                    settings['bar_border_color'] = st.color_picker("Bar Border",
                                                                   str(settings.get('bar_border_color', 'black')),
                                                                   key=f"bordc_{gene}")
                    settings['error_color'] = st.color_picker("Error Bars",
                                                             str(settings.get('error_color', 'black')),
                                                             key=f"errc_{gene}")
                
                with col3:
                    settings['plot_bgcolor'] = st.color_picker("Plot Background",
                                                               str(settings.get('plot_bgcolor', 'white')),
                                                               key=f"plotbg_{gene}")
                    settings['reference_line_color'] = st.color_picker("Reference Line",
                                                                       str(settings.get('reference_line_color', 'gray')),
                                                                       key=f"refc_{gene}")
            
            # Generate and display graph
            st.markdown("### 📊 Live Preview")
            
            fig = GraphGenerator.create_graph(gene_data, gene, settings)
            st.plotly_chart(fig, width='stretch', key=f"graph_{gene}")
            
            if 'graphs' not in st.session_state:
                st.session_state.graphs = {}
            st.session_state.graphs[gene] = fig
            
            st.markdown("---")

# ==================== TAB 5: EXPORT ====================
with tab5:
    st.header("Export Results")
    
    if st.session_state.processed_data:
        params = {
            'Date': datetime.now().strftime("%Y-%m-%d %H:%M"),
            'Efficacy': st.session_state.selected_efficacy,
            'HK_Gene': st.session_state.hk_gene,
            'Baseline': baseline_condition if 'baseline_condition' in locals() else 'N/A',
            'Pvalue_Ref': pvalue_condition if 'pvalue_condition' in locals() else 'N/A',
            'Genes': len(st.session_state.processed_data)
        }
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("### 📊 Excel Report")
            excel = export_excel(st.session_state.data, st.session_state.processed_data,
                                params, st.session_state.sample_mapping)
            
            st.download_button("📥 Download Excel",
                             excel,
                             f"qPCR_{st.session_state.selected_efficacy}_{datetime.now().strftime('%Y%m%d')}.xlsx",
                             "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                             type="primary")
        
        with col2:
            st.markdown("### 📈 All Graphs (HTML)")
            
            if st.session_state.graphs:
                html = ["<html><head><title>qPCR Graphs</title></head><body>",
                       f"<h1>{st.session_state.selected_efficacy} Analysis</h1>"]
                
                for gene, fig in st.session_state.graphs.items():
                    html.append(f"<h2>{gene}</h2>")
                    html.append(fig.to_html(include_plotlyjs='cdn'))
                    html.append("<hr>")
                
                html.append("</body></html>")
                
                st.download_button("📥 Download Graphs",
                                 "\n".join(html),
                                 f"graphs_{datetime.now().strftime('%Y%m%d')}.html",
                                 "text/html",
                                 type="primary")
        
        st.markdown("---")
        st.subheader("Individual Files")
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.markdown("**CSV per Gene**")
            for gene, df in st.session_state.processed_data.items():
                csv = io.StringIO()
                df.drop(columns=['Target_Ct_Values'], errors='ignore').to_csv(csv, index=False)
                st.download_button(f"📥 {gene}.csv",
                                 csv.getvalue(),
                                 f"{gene}.csv",
                                 "text/csv",
                                 key=f"csv_{gene}")
        
        with col2:
            st.markdown("**HTML per Gene**")
            for gene, fig in st.session_state.graphs.items():
                html = io.StringIO()
                fig.write_html(html)
                st.download_button(f"📥 {gene}.html",
                                 html.getvalue(),
                                 f"{gene}.html",
                                 "text/html",
                                 key=f"html_{gene}")
        
        with col3:
            st.markdown("**Config**")
            config = {'params': params, 'mapping': st.session_state.sample_mapping,
                     'settings': st.session_state.graph_settings}
            st.download_button("📥 Config.json",
                             json.dumps(config, indent=2),
                             f"config_{datetime.now().strftime('%Y%m%d')}.json",
                             "application/json")

st.markdown("---")
st.markdown("<div style='text-align: center; color: #666;'><p>🧬 qPCR Analysis Pro v4.0 | Optimized & Production-Ready</p></div>", 
           unsafe_allow_html=True)
