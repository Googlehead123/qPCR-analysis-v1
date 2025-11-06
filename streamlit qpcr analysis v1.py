import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from scipy import stats
import io
import json
from datetime import datetime
from typing import Dict, List, Tuple, Optional

# ==================== PAGE CONFIG ====================
st.set_page_config(page_title="qPCR Analysis Suite Ultimate", layout="wide", initial_sidebar_state="expanded")

# ==================== SESSION STATE INIT ====================
defaults = {
    'data': None, 'processed_data': {}, 'sample_mapping': {}, 'analysis_templates': {},
    'graphs': {}, 'excluded_wells': set(), 'excluded_samples': set(), 'selected_efficacy': None,
    'hk_gene': None, 'included_samples_for_analysis': {}, 'graph_settings': {}
}
for key, val in defaults.items():
    if key not in st.session_state:
        st.session_state[key] = val

# ==================== EFFICACY DATABASE ====================
EFFICACY_CONFIG = {
    '탄력': {
        'genes': ['COL1A1', 'ELN', 'FBN-1', 'FBN1'], 'cell': 'HS68 fibroblast',
        'controls': {'baseline': 'Non-treated', 'positive': 'TGFb'},
        'description': 'Elasticity - All vs Non-treated baseline'
    },
    '항노화': {
        'genes': ['COL1A1', 'COL1', 'MMP-1', 'MMP1'], 'cell': 'HS68 fibroblast',
        'controls': {'baseline': 'Non-treated (No UV)', 'negative': 'UVB only', 'positive': 'UVB+TGFb'},
        'description': 'Anti-aging - COL1↑ MMP1↓ vs Non-treated baseline',
        'expected_direction': {'COL1A1': 'up', 'COL1': 'up', 'MMP-1': 'down', 'MMP1': 'down'}
    },
    '보습': {
        'genes': ['AQP3', 'HAS3'], 'cell': 'HaCaT keratinocyte',
        'controls': {'baseline': 'Non-treated', 'positive': 'Retinoic acid'},
        'description': 'Hydration - All vs Non-treated baseline'
    },
    '장벽': {
        'genes': ['FLG', 'CLDN', 'IVL'], 'cell': 'HaCaT keratinocyte',
        'controls': {'baseline': 'Non-treated', 'positive': 'Retinoic acid'},
        'description': 'Barrier - All vs Non-treated baseline'
    },
    '표피증식': {
        'genes': ['KI67', 'PCNA'], 'cell': 'HaCaT keratinocyte',
        'controls': {'baseline': 'Non-treated', 'positive': 'TGFb or FBS'},
        'description': 'Proliferation - All vs Non-treated baseline'
    },
    '멜라닌억제': {
        'genes': ['MITF', 'TYR', 'Melanin'], 'cell': 'B16F10 melanocyte',
        'controls': {'baseline': 'Non-treated', 'negative': 'α-MSH only', 'positive': 'α-MSH+Arbutin'},
        'description': 'Melanin inhibition - All vs Non-treated baseline',
        'expected_direction': {'MITF': 'down', 'TYR': 'down', 'Melanin': 'down'}
    },
    '진정': {
        'genes': ['IL1B', 'IL-1β', 'IL6', 'TNFA', 'TNFα'], 'cell': 'HaCaT keratinocyte',
        'controls': {'baseline': 'Non-treated', 'negative': 'IL4+PolyIC', 'positive': 'Inflammation+Dex'},
        'description': 'Anti-inflammation - All vs Non-treated baseline',
        'expected_direction': {'IL1B': 'down', 'IL-1β': 'down', 'IL6': 'down', 'TNFA': 'down', 'TNFα': 'down'}
    },
    '지질억제': {
        'genes': ['SREBPA', 'SREBPa', 'SREBPC', 'SREBPc', 'PPARY', 'PPARy'], 'cell': 'SZ95 sebocyte',
        'controls': {'baseline': 'Non-treated', 'negative': 'IGF only', 'positive': 'IGF+inhibitor'},
        'description': 'Sebum inhibition - All vs Non-treated baseline',
        'expected_direction': {'SREBPA': 'down', 'SREBPa': 'down', 'SREBPC': 'down', 'SREBPc': 'down'}
    },
    '냉감': {
        'genes': ['TRPM8', 'CIRBP'], 'cell': 'HaCaT keratinocyte',
        'controls': {'baseline': 'Non-treated', 'positive': 'Menthol'},
        'description': 'Cooling - All vs Non-treated baseline'
    }
}

# ==================== PARSER CLASS ====================
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
    def parse_format1(df, start):
        df = df.iloc[start:].reset_index(drop=True)
        df.columns = df.iloc[0]
        df = df.iloc[1:].reset_index(drop=True)
        
        well_col = next((c for c in ['Well Position', 'Well'] if c in df.columns), df.columns[0])
        ct_col = next((c for c in ['CT', 'Ct', 'Cт'] if c in df.columns), None)
        
        if not ct_col:
            return None
        
        parsed = pd.DataFrame({
            'Well': df[well_col],
            'Sample': df.get('Sample Name', df.iloc[:, 2]),
            'Target': df.get('Target Name', df.iloc[:, 3]),
            'CT': pd.to_numeric(df[ct_col], errors='coerce')
        })
        
        return parsed.dropna(subset=['CT']).query('Sample.notna() & Target.notna()')
    
    @staticmethod
    def parse_format2(df, start):
        df = df.iloc[start:].reset_index(drop=True)
        df.columns = df.iloc[0]
        df = df.iloc[1:].reset_index(drop=True)
        
        return pd.DataFrame({
            'Well': df['Well'],
            'Sample': df['Sample Name'],
            'Target': df['Target Name'],
            'CT': pd.to_numeric(df['Cт'], errors='coerce')
        }).dropna(subset=['CT']).query('Sample.notna() & Target.notna()')
    
    @staticmethod
    def parse(file):
        try:
            df = None
            for enc in ['utf-8', 'latin-1', 'cp1252']:
                try:
                    df = pd.read_csv(file, encoding=enc, low_memory=False, skip_blank_lines=False)
                    break
                except UnicodeDecodeError:
                    continue
            
            if df is None:
                return None
            
            fmt, start = QPCRParser.detect_format(df)
            return QPCRParser.parse_format1(df, start) if fmt == 'format1' else QPCRParser.parse_format2(df, start) if fmt == 'format2' else None
        except Exception as e:
            st.error(f"Parse error: {e}")
            return None

# ==================== ANALYSIS ENGINE ====================
class AnalysisEngine:
    @staticmethod
    def calculate_ddct(data: pd.DataFrame, hk_gene: str, baseline_condition: str, 
                       excluded_wells: set, excluded_samples: set, sample_mapping: dict,
                       included_samples_per_gene: Dict[str, List[str]]) -> Dict[str, pd.DataFrame]:
        """Gene-by-gene ΔΔCt with baseline as reference"""
        
        # Filter data
        data = data[~data['Well'].isin(excluded_wells) & ~data['Sample'].isin(excluded_samples)].copy()
        
        # Apply sample name mapping
        data['Condition'] = data['Sample'].map(lambda x: sample_mapping.get(x, {}).get('condition', x))
        data['Group'] = data['Sample'].map(lambda x: sample_mapping.get(x, {}).get('group', 'Treatment'))
        
        gene_results = {}
        
        # Process each target gene separately
        for target in data['Target'].unique():
            if target.upper() in [hk_gene.upper(), 'ACTIN', 'B-ACTIN', 'GAPDH', 'ACTB']:
                continue
            
            # Filter by included samples for this gene
            included_samples = included_samples_per_gene.get(target, [])
            if not included_samples:
                continue
            
            target_data = data[data['Target'] == target]
            target_data = target_data[target_data['Condition'].isin(included_samples)]
            
            # Get baseline reference
            baseline_target = target_data[target_data['Condition'] == baseline_condition]
            baseline_hk = data[(data['Condition'] == baseline_condition) & (data['Target'] == hk_gene)]
            
            if len(baseline_target) == 0 or len(baseline_hk) == 0:
                continue
            
            baseline_delta_ct = baseline_target['CT'].mean() - baseline_hk['CT'].mean()
            
            results = []
            
            for condition in target_data['Condition'].unique():
                cond_data = target_data[target_data['Condition'] == condition]
                hk_data = data[(data['Condition'] == condition) & (data['Target'] == hk_gene)]
                
                if len(hk_data) == 0:
                    continue
                
                # Get raw CT values for t-test
                target_ct_values = cond_data['CT'].values
                hk_ct_values = hk_data['CT'].values
                
                # Calculate ΔCt = Target_Ct - HK_Ct
                target_ct_mean = target_ct_values.mean()
                hk_ct_mean = hk_ct_values.mean()
                delta_ct = target_ct_mean - hk_ct_mean
                
                # ΔΔCt relative to baseline
                ddct = delta_ct - baseline_delta_ct
                rel_expr = 2 ** (-ddct)
                
                # Calculate SEM
                ct_sd = target_ct_values.std() if len(target_ct_values) > 1 else 0
                sem = ct_sd / np.sqrt(len(target_ct_values)) if len(target_ct_values) > 1 else 0
                
                # Store raw values for t-test
                original_sample = cond_data['Sample'].iloc[0]
                group = sample_mapping.get(original_sample, {}).get('group', 'Treatment')
                
                results.append({
                    'Target': target,
                    'Condition': condition,
                    'Original_Sample': original_sample,
                    'Group': group,
                    'n_replicates': len(target_ct_values),
                    'Target_Ct_Mean': target_ct_mean,
                    'Target_Ct_SD': ct_sd,
                    'Target_Ct_Values': target_ct_values,  # Store for t-test
                    'HK_Ct_Mean': hk_ct_mean,
                    'Delta_Ct': delta_ct,
                    'Delta_Delta_Ct': ddct,
                    'Relative_Expression': rel_expr,
                    'SEM': sem,
                    'Fold_Change': rel_expr
                })
            
            if results:
                gene_results[target] = pd.DataFrame(results)
        
        return gene_results
    
    @staticmethod
    def calculate_statistics(gene_data: pd.DataFrame, reference_condition: str) -> pd.DataFrame:
        """Two-tailed t-test vs reference (Excel TTEST equivalent)"""
        results = gene_data.copy()
        results['p_value'] = np.nan
        results['significance'] = ''
        
        # Get reference values
        ref_data = results[results['Condition'] == reference_condition]
        if len(ref_data) == 0:
            return results
        
        ref_values = ref_data.iloc[0]['Target_Ct_Values']
        
        for idx, row in results.iterrows():
            if row['Condition'] == reference_condition:
                results.at[idx, 'p_value'] = 1.0  # Same as reference
                continue
            
            sample_values = row['Target_Ct_Values']
            
            # Two-tailed independent t-test (type 2, 2 in Excel TTEST)
            try:
                if len(ref_values) > 1 and len(sample_values) > 1:
                    t_stat, p_val = stats.ttest_ind(ref_values, sample_values)
                elif len(ref_values) == 1 and len(sample_values) > 1:
                    t_stat, p_val = stats.ttest_1samp(sample_values, ref_values[0])
                elif len(sample_values) == 1 and len(ref_values) > 1:
                    t_stat, p_val = stats.ttest_1samp(ref_values, sample_values[0])
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
    def create_gene_graph(data: pd.DataFrame, gene: str, settings: dict, 
                         included_conditions: List[str], efficacy_config: dict = None) -> go.Figure:
        """Highly customizable gene graph"""
        
        # Filter by included conditions
        gene_data = data[data['Condition'].isin(included_conditions)].copy()
        
        # Custom sorting if specified
        if settings.get('custom_order'):
            order_map = {cond: i for i, cond in enumerate(settings['custom_order'])}
            gene_data['sort_key'] = gene_data['Condition'].map(lambda x: order_map.get(x, 999))
            gene_data = gene_data.sort_values('sort_key')
        
        fig = go.Figure()
        
        # Get colors for each condition
        colors = [settings['condition_colors'].get(cond, '#636EFA') for cond in gene_data['Condition']]
        
        # Create bar chart
        fig.add_trace(go.Bar(
            x=gene_data['Condition'] if settings.get('orientation', 'vertical') == 'vertical' else gene_data['Fold_Change'],
            y=gene_data['Fold_Change'] if settings.get('orientation', 'vertical') == 'vertical' else gene_data['Condition'],
            orientation='v' if settings.get('orientation', 'vertical') == 'vertical' else 'h',
            error_y=dict(
                type='data',
                array=gene_data['SEM'] * 1.96 if settings.get('show_error', True) else None,
                visible=settings.get('show_error', True),
                thickness=settings.get('error_bar_thickness', 2),
                width=settings.get('error_bar_width', 4)
            ) if settings.get('orientation', 'vertical') == 'vertical' else None,
            error_x=dict(
                type='data',
                array=gene_data['SEM'] * 1.96 if settings.get('show_error', True) else None,
                visible=settings.get('show_error', True),
                thickness=settings.get('error_bar_thickness', 2),
                width=settings.get('error_bar_width', 4)
            ) if settings.get('orientation', 'vertical') == 'horizontal' else None,
            text=gene_data['significance'] if settings.get('show_significance', True) else None,
            textposition=settings.get('sig_position', 'outside'),
            textfont=dict(size=settings.get('sig_font_size', 16), color=settings.get('sig_color', 'black')),
            marker=dict(
                color=colors,
                line=dict(color=settings.get('bar_border_color', 'black'), 
                         width=settings.get('bar_border_width', 0))
            ),
            width=settings.get('bar_width', 0.8),
            showlegend=False
        ))
        
        # Add reference line
        if settings.get('show_reference_line', True):
            ref_val = settings.get('reference_value', 1.0)
            if settings.get('orientation', 'vertical') == 'vertical':
                fig.add_hline(y=ref_val, line_dash="dash", line_color="gray", 
                             annotation_text=f"Reference ({ref_val})",
                             annotation_position="right")
            else:
                fig.add_vline(x=ref_val, line_dash="dash", line_color="gray",
                             annotation_text=f"Reference ({ref_val})")
        
        # Add expected direction annotation
        if efficacy_config and 'expected_direction' in efficacy_config:
            direction = efficacy_config['expected_direction'].get(gene, '')
            if direction and settings.get('show_expected_direction', True):
                direction_text = '↑ Expected increase' if direction == 'up' else '↓ Expected decrease'
                fig.add_annotation(
                    text=direction_text,
                    xref='paper', yref='paper',
                    x=settings.get('direction_annotation_x', 0.02),
                    y=settings.get('direction_annotation_y', 0.98),
                    showarrow=False,
                    font=dict(size=settings.get('direction_font_size', 12), 
                             color=settings.get('direction_color', 'red')),
                    align='left'
                )
        
        # Layout customization
        fig.update_layout(
            title=dict(
                text=settings.get('title', f"{gene} Expression"),
                font=dict(size=settings.get('title_size', 20), color=settings.get('title_color', 'black')),
                x=settings.get('title_position_x', 0.5),
                y=settings.get('title_position_y', 0.95),
                xanchor='center'
            ),
            xaxis=dict(
                title=dict(text=settings.get('xlabel', 'Condition'),
                          font=dict(size=settings.get('axis_label_size', 14), 
                                   color=settings.get('axis_label_color', 'black'))),
                showgrid=settings.get('show_grid', True),
                gridcolor=settings.get('grid_color', 'lightgray'),
                tickangle=settings.get('x_tick_angle', 0),
                tickfont=dict(size=settings.get('tick_font_size', 12))
            ),
            yaxis=dict(
                title=dict(text=settings.get('ylabel', 'Fold Change'),
                          font=dict(size=settings.get('axis_label_size', 14),
                                   color=settings.get('axis_label_color', 'black'))),
                showgrid=settings.get('show_grid', True),
                gridcolor=settings.get('grid_color', 'lightgray'),
                tickfont=dict(size=settings.get('tick_font_size', 12)),
                range=settings.get('y_range', None)
            ),
            template=settings.get('color_scheme', 'plotly_white'),
            height=settings.get('figure_height', 600),
            width=settings.get('figure_width', 1000),
            plot_bgcolor=settings.get('plot_bgcolor', 'white'),
            paper_bgcolor=settings.get('paper_bgcolor', 'white'),
            margin=dict(
                l=settings.get('margin_left', 80),
                r=settings.get('margin_right', 80),
                t=settings.get('margin_top', 100),
                b=settings.get('margin_bottom', 80)
            )
        )
        
        return fig

# ==================== EXPORT ====================
def export_to_excel(raw_data: pd.DataFrame, processed_data: Dict[str, pd.DataFrame], 
                   params: dict, mapping: dict) -> bytes:
    output = io.BytesIO()
    
    with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
        pd.DataFrame([params]).to_excel(writer, sheet_name='Parameters', index=False)
        pd.DataFrame([{'Original': k, **v} for k, v in mapping.items()]).to_excel(
            writer, sheet_name='Sample_Mapping', index=False
        )
        raw_data.to_excel(writer, sheet_name='Raw_Data', index=False)
        
        for gene, gene_df in processed_data.items():
            # Remove raw CT values column before export
            export_df = gene_df.drop(columns=['Target_Ct_Values'], errors='ignore')
            export_df.to_excel(writer, sheet_name=f"{gene}"[:31], index=False)
        
        if processed_data:
            all_data = pd.concat([df.drop(columns=['Target_Ct_Values'], errors='ignore') 
                                 for df in processed_data.values()], ignore_index=True)
            summary = all_data.groupby(['Target', 'Group']).agg({
                'Relative_Expression': ['mean', 'std', 'count'],
                'p_value': 'min'
            }).round(4)
            summary.to_excel(writer, sheet_name='Summary')
    
    return output.getvalue()

# ==================== UI ====================
st.title("🧬 qPCR Analysis Suite Ultimate")
st.markdown("**Complete control: Filter • Analyze • Customize • Export**")

# Sidebar
with st.sidebar:
    st.header("📖 Quick Guide")
    st.markdown("""
    **Workflow:**
    1. Upload CSV
    2. Filter samples & genes
    3. Select efficacy type
    4. Map samples
    5. Choose baseline reference
    6. Analyze gene-by-gene
    7. Select samples per graph
    8. Customize graphs fully
    9. Export everything
    """)
    
    st.subheader("📋 Templates")
    template_name = st.text_input("Template name:")
    if st.button("💾 Save") and st.session_state.sample_mapping:
        st.session_state.analysis_templates[template_name] = {
            'mapping': st.session_state.sample_mapping.copy(),
            'efficacy': st.session_state.selected_efficacy,
            'graph_settings': st.session_state.graph_settings.copy(),
            'timestamp': datetime.now().isoformat()
        }
        st.success(f"✅ '{template_name}' saved")
    
    if st.session_state.analysis_templates:
        load = st.selectbox("Load:", [""] + list(st.session_state.analysis_templates.keys()))
        if load:
            t = st.session_state.analysis_templates[load]
            st.session_state.sample_mapping = t['mapping']
            st.session_state.selected_efficacy = t.get('efficacy')
            st.session_state.graph_settings = t.get('graph_settings', {})
            st.info(f"✅ Loaded '{load}'")

# Tabs
tab1, tab2, tab3, tab4, tab5 = st.tabs(["📁 Data", "🗺️ Mapping", "🔬 Analysis", "📊 Graphs", "📤 Export"])

# ==================== TAB 1: DATA ====================
with tab1:
    st.header("Step 1: Upload & Filter Data")
    
    uploaded_files = st.file_uploader("Upload qPCR CSV", type=['csv'], accept_multiple_files=True)
    
    if uploaded_files:
        all_data = []
        for file in uploaded_files:
            parsed = QPCRParser.parse(file)
            if parsed is not None:
                parsed['Source_File'] = file.name
                all_data.append(parsed)
                st.success(f"✅ {file.name}: {len(parsed)} wells")
        
        if all_data:
            st.session_state.data = pd.concat(all_data, ignore_index=True)
            
            col1, col2, col3, col4 = st.columns(4)
            col1.metric("Wells", len(st.session_state.data))
            col2.metric("Samples", st.session_state.data['Sample'].nunique())
            col3.metric("Genes", st.session_state.data['Target'].nunique())
            
            # Detect HK gene
            hk_genes = [g for g in st.session_state.data['Target'].unique() 
                       if g.upper() in ['ACTIN', 'B-ACTIN', 'GAPDH', 'ACTB']]
            if hk_genes:
                st.session_state.hk_gene = st.selectbox("🔬 Housekeeping Gene", hk_genes)
                col4.metric("HK Gene", st.session_state.hk_gene)
            
            # Sample filter
            st.subheader("🎯 Sample Filter")
            all_samples = sorted(st.session_state.data['Sample'].unique())
            selected_samples = st.multiselect(
                "Include samples:",
                all_samples,
                default=all_samples
            )
            st.session_state.excluded_samples = set(all_samples) - set(selected_samples)
            
            # Gene filter
            st.subheader("🧬 Gene Filter")
            all_genes = [g for g in sorted(st.session_state.data['Target'].unique()) 
                        if g.upper() not in ['ACTIN', 'B-ACTIN', 'GAPDH', 'ACTB']]
            selected_genes = st.multiselect(
                "Include genes:",
                all_genes,
                default=all_genes
            )
            
            # Preview
            display_data = st.session_state.data[
                st.session_state.data['Sample'].isin(selected_samples) &
                st.session_state.data['Target'].isin(selected_genes + [st.session_state.hk_gene])
            ]
            
            st.subheader("📊 Filtered Preview")
            st.dataframe(display_data.head(50), height=300)
            
            # Outliers
            st.subheader("🚩 Outlier Detection")
            for target in selected_genes:
                target_data = display_data[display_data['Target'] == target]
                q1, q3 = target_data['CT'].quantile([0.25, 0.75])
                iqr = q3 - q1
                outliers = target_data[(target_data['CT'] < q1 - 1.5*iqr) | (target_data['CT'] > q3 + 1.5*iqr)]
                
                if len(outliers) > 0:
                    with st.expander(f"⚠️ {target}: {len(outliers)} outliers"):
                        st.dataframe(outliers[['Well', 'Sample', 'CT']])
                        if st.checkbox(f"Exclude", key=f"out_{target}"):
                            st.session_state.excluded_wells.update(outliers['Well'].tolist())

# ==================== TAB 2: MAPPING ====================
with tab2:
    st.header("Step 2: Sample Mapping")
    
    if st.session_state.data is not None:
        # Efficacy selection
        detected = set(st.session_state.data['Target'].unique())
        suggested = next((e for e, c in EFFICACY_CONFIG.items() if any(g in detected for g in c['genes'])), None)
        
        efficacy = st.selectbox(
            "🎯 Efficacy Type",
            list(EFFICACY_CONFIG.keys()),
            index=list(EFFICACY_CONFIG.keys()).index(suggested) if suggested else 0
        )
        st.session_state.selected_efficacy = efficacy
        
        config = EFFICACY_CONFIG[efficacy]
        st.info(f"**{config['description']}**")
        
        with st.expander("📋 Control Structure"):
            for k, v in config['controls'].items():
                st.markdown(f"- **{k.title()}**: {v}")
        
        # Mapping interface
        st.subheader("🗺️ Condition Mapping")
        
        samples = [s for s in sorted(st.session_state.data['Sample'].unique()) 
                  if s not in st.session_state.excluded_samples]
        
        group_types = ['Baseline', 'Negative Control', 'Positive Control', 'Treatment']
        
        for sample in samples:
            if sample not in st.session_state.sample_mapping:
                st.session_state.sample_mapping[sample] = {
                    'condition': sample, 'group': 'Treatment', 'concentration': ''
                }
            
            col1, col2, col3, col4 = st.columns([1, 3, 2, 2])
            
            with col1:
                st.text(sample)
            with col2:
                cond = st.text_input("Condition", st.session_state.sample_mapping[sample]['condition'],
                                    key=f"c_{sample}", label_visibility="collapsed")
                st.session_state.sample_mapping[sample]['condition'] = cond
            with col3:
                grp = st.selectbox("Group", group_types,
                                  index=group_types.index(st.session_state.sample_mapping[sample]['group']) 
                                  if st.session_state.sample_mapping[sample]['group'] in group_types else 0,
                                  key=f"g_{sample}", label_visibility="collapsed")
                st.session_state.sample_mapping[sample]['group'] = grp
            with col4:
                conc = st.text_input("Conc.", st.session_state.sample_mapping[sample]['concentration'],
                                    key=f"conc_{sample}", label_visibility="collapsed")
                st.session_state.sample_mapping[sample]['concentration'] = conc
        
        st.subheader("📊 Mapping Summary")
        mapping_df = pd.DataFrame([{'Original': k, **v} for k, v in st.session_state.sample_mapping.items()])
        st.dataframe(mapping_df, use_container_width=True)
    else:
        st.warning("⚠️ Upload data first")

# ==================== TAB 3: ANALYSIS ====================
with tab3:
    st.header("Step 3: Analysis Setup & Execution")
    
    if st.session_state.data is not None and st.session_state.sample_mapping and st.session_state.hk_gene:
        config = EFFICACY_CONFIG.get(st.session_state.selected_efficacy, {})
        
        # Select baseline reference
        st.subheader("📌 Reference Selection")
        
        col1, col2 = st.columns(2)
        
        with col1:
            baseline_samples = [k for k, v in st.session_state.sample_mapping.items() 
                              if v['group'] == 'Baseline']
            if not baseline_samples:
                baseline_samples = [k for k, v in st.session_state.sample_mapping.items() 
                                  if 'non-treated' in v['condition'].lower() or 'control' in v['condition'].lower()]
            
            baseline_sample = st.selectbox(
                "🎯 Baseline Reference (Fold Change = 1.0)",
                baseline_samples if baseline_samples else list(st.session_state.sample_mapping.keys()),
                help="All fold changes calculated relative to this baseline"
            )
            baseline_condition = st.session_state.sample_mapping[baseline_sample]['condition']
        
        with col2:
            pvalue_samples = [k for k, v in st.session_state.sample_mapping.items()]
            pvalue_sample = st.selectbox(
                "📊 P-value Reference (Statistical Comparison)",
                pvalue_samples,
                index=pvalue_samples.index(baseline_sample) if baseline_sample in pvalue_samples else 0,
                help="P-values calculated comparing each sample to this reference"
            )
            pvalue_condition = st.session_state.sample_mapping[pvalue_sample]['condition']
        
        st.info(f"**Setup:** Fold changes vs **{baseline_condition}** | P-values vs **{pvalue_condition}**")
        
        # Per-gene sample inclusion
        st.subheader("🧬 Select Samples Per Gene")
        st.markdown("*Choose which samples to include in analysis for each gene*")
        
        # Get available conditions
        available_conditions = [v['condition'] for v in st.session_state.sample_mapping.values()]
        
        # Get genes from filtered data
        genes = [g for g in st.session_state.data['Target'].unique() 
                if g.upper() not in ['ACTIN', 'B-ACTIN', 'GAPDH', 'ACTB']]
        
        if 'included_samples_for_analysis' not in st.session_state:
            st.session_state.included_samples_for_analysis = {}
        
        for gene in genes:
            if gene not in st.session_state.included_samples_for_analysis:
                st.session_state.included_samples_for_analysis[gene] = available_conditions.copy()
            
            with st.expander(f"📍 {gene} - Sample Selection", expanded=False):
                selected = st.multiselect(
                    f"Include these conditions for {gene}:",
                    available_conditions,
                    default=st.session_state.included_samples_for_analysis[gene],
                    key=f"incl_{gene}"
                )
                st.session_state.included_samples_for_analysis[gene] = selected
                
                if baseline_condition not in selected:
                    st.warning(f"⚠️ Baseline '{baseline_condition}' not included! Add it for proper analysis.")
        
        # Run analysis
        if st.button("🔬 Run Analysis", type="primary"):
            with st.spinner("Calculating ΔΔCt for each gene..."):
                gene_results = AnalysisEngine.calculate_ddct(
                    st.session_state.data,
                    st.session_state.hk_gene,
                    baseline_condition,
                    st.session_state.excluded_wells,
                    st.session_state.excluded_samples,
                    st.session_state.sample_mapping,
                    st.session_state.included_samples_for_analysis
                )
                
                if gene_results:
                    # Calculate statistics for each gene
                    for gene, gene_df in gene_results.items():
                        gene_results[gene] = AnalysisEngine.calculate_statistics(gene_df, pvalue_condition)
                    
                    st.session_state.processed_data = gene_results
                    st.success(f"✅ {len(gene_results)} genes analyzed!")
        
        # Display results
        if st.session_state.processed_data:
            st.subheader("📊 Analysis Results")
            
            all_results = pd.concat([df.drop(columns=['Target_Ct_Values'], errors='ignore') 
                                    for df in st.session_state.processed_data.values()], ignore_index=True)
            
            col1, col2, col3 = st.columns(3)
            col1.metric("Genes", len(st.session_state.processed_data))
            col2.metric("Conditions", all_results['Condition'].nunique())
            sig_count = (all_results['p_value'] < 0.05).sum()
            col3.metric("Significant (p<0.05)", f"{sig_count}/{len(all_results)}")
            
            for gene, gene_df in st.session_state.processed_data.items():
                with st.expander(f"🧬 {gene} Results", expanded=False):
                    config = EFFICACY_CONFIG.get(st.session_state.selected_efficacy, {})
                    if 'expected_direction' in config and gene in config['expected_direction']:
                        direction = config['expected_direction'][gene]
                        st.caption(f"Expected: {'↑ Increase' if direction == 'up' else '↓ Decrease'}")
                    
                    display_cols = ['Condition', 'Group', 'Fold_Change', 'p_value', 'significance', 
                                  'n_replicates', 'Target_Ct_Mean', 'HK_Ct_Mean', 'Delta_Ct', 'SEM']
                    
                    display_df = gene_df[[c for c in display_cols if c in gene_df.columns]]
                    
                    styled = display_df.style.background_gradient(
                        subset=['Fold_Change'], cmap='RdYlGn', vmin=0, vmax=3
                    ).format({
                        'Fold_Change': '{:.3f}',
                        'p_value': '{:.4f}',
                        'Target_Ct_Mean': '{:.2f}',
                        'HK_Ct_Mean': '{:.2f}',
                        'Delta_Ct': '{:.2f}',
                        'SEM': '{:.3f}'
                    }, na_rep='—')
                    
                    st.dataframe(styled, use_container_width=True)
    else:
        st.warning("⚠️ Complete previous steps")

# ==================== TAB 4: GRAPHS ====================
with tab4:
    st.header("Step 4: Interactive Graph Customization")
    
    if st.session_state.processed_data:
        # Initialize graph settings per gene
        for gene in st.session_state.processed_data.keys():
            if gene not in st.session_state.graph_settings:
                default_colors = px.colors.qualitative.Plotly
                idx = list(st.session_state.processed_data.keys()).index(gene)
                
                st.session_state.graph_settings[gene] = {
                    'title': f"{gene} Expression",
                    'xlabel': 'Condition',
                    'ylabel': 'Fold Change (Relative to Baseline)',
                    'title_size': 20,
                    'font_size': 14,
                    'sig_font_size': 16,
                    'figure_width': 1000,
                    'figure_height': 600,
                    'color_scheme': 'plotly_white',
                    'show_error': True,
                    'show_significance': True,
                    'show_grid': True,
                    'orientation': 'vertical',
                    'bar_width': 0.8,
                    'show_reference_line': True,
                    'reference_value': 1.0,
                    'show_expected_direction': True,
                    'condition_colors': {},
                    'included_conditions': list(st.session_state.processed_data[gene]['Condition'].unique()),
                    'custom_order': None,
                    'x_tick_angle': -45,
                    'tick_font_size': 12,
                    'axis_label_size': 14,
                    'title_color': 'black',
                    'axis_label_color': 'black',
                    'sig_color': 'black',
                    'sig_position': 'outside',
                    'bar_border_width': 0,
                    'bar_border_color': 'black',
                    'error_bar_thickness': 2,
                    'error_bar_width': 4,
                    'grid_color': 'lightgray',
                    'plot_bgcolor': 'white',
                    'paper_bgcolor': 'white',
                    'y_range': None,
                    'margin_left': 80,
                    'margin_right': 80,
                    'margin_top': 100,
                    'margin_bottom': 80,
                    'title_position_x': 0.5,
                    'title_position_y': 0.95,
                    'direction_annotation_x': 0.02,
                    'direction_annotation_y': 0.98,
                    'direction_font_size': 12,
                    'direction_color': 'red'
                }
        
        # Global presets
        st.subheader("⚡ Quick Presets (Apply to All Genes)")
        col1, col2, col3, col4 = st.columns(4)
        
        if col1.button("📄 Publication Style"):
            for gene in st.session_state.graph_settings:
                st.session_state.graph_settings[gene].update({
                    'color_scheme': 'simple_white', 'font_size': 14, 'title_size': 18,
                    'figure_width': 1000, 'figure_height': 600, 'bar_border_width': 1
                })
            st.rerun()
        
        if col2.button("🎨 Presentation Style"):
            for gene in st.session_state.graph_settings:
                st.session_state.graph_settings[gene].update({
                    'color_scheme': 'presentation', 'font_size': 18, 'title_size': 24,
                    'figure_width': 1200, 'figure_height': 700
                })
            st.rerun()
        
        if col3.button("🌙 Dark Mode"):
            for gene in st.session_state.graph_settings:
                st.session_state.graph_settings[gene].update({
                    'color_scheme': 'plotly_dark', 'plot_bgcolor': '#111', 'paper_bgcolor': '#111',
                    'title_color': 'white', 'axis_label_color': 'white', 'sig_color': 'white'
                })
            st.rerun()
        
        if col4.button("🔄 Reset All"):
            st.session_state.graph_settings = {}
            st.rerun()
        
        st.markdown("---")
        
        # Gene-specific customization
        efficacy_config = EFFICACY_CONFIG.get(st.session_state.selected_efficacy, {})
        
        for gene in st.session_state.processed_data.keys():
            st.markdown(f"## 🧬 {gene} Graph")
            
            settings = st.session_state.graph_settings[gene]
            gene_data = st.session_state.processed_data[gene]
            
            # Customization tabs
            tab_samples, tab_style, tab_colors, tab_labels, tab_advanced = st.tabs([
                "📋 Samples", "🎨 Style", "🌈 Colors", "📝 Labels", "⚙️ Advanced"
            ])
            
            # TAB: Sample Selection
            with tab_samples:
                st.markdown("**Select samples to display in this graph**")
                
                all_conditions = list(gene_data['Condition'].unique())
                settings['included_conditions'] = st.multiselect(
                    "Include conditions:",
                    all_conditions,
                    default=settings['included_conditions'],
                    key=f"samples_{gene}"
                )
                
                st.markdown("**Custom order (optional)**")
                if st.checkbox(f"Use custom order for {gene}", key=f"order_{gene}"):
                    settings['custom_order'] = st.multiselect(
                        "Drag to reorder:",
                        settings['included_conditions'],
                        default=settings['included_conditions'],
                        key=f"custom_order_{gene}"
                    )
                else:
                    settings['custom_order'] = None
            
            # TAB: Style
            with tab_style:
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    settings['orientation'] = st.radio("Orientation", ['vertical', 'horizontal'], 
                                                      key=f"orient_{gene}")
                    settings['color_scheme'] = st.selectbox("Theme", 
                        ['plotly_white', 'plotly', 'plotly_dark', 'seaborn', 'simple_white', 
                         'presentation', 'ggplot2', 'none'],
                        key=f"theme_{gene}")
                    settings['bar_width'] = st.slider("Bar Width", 0.1, 1.5, settings['bar_width'], 
                                                     key=f"barw_{gene}")
                
                with col2:
                    settings['show_error'] = st.checkbox("Error Bars", settings['show_error'], 
                                                        key=f"err_{gene}")
                    if settings['show_error']:
                        settings['error_bar_thickness'] = st.slider("Error Bar Thickness", 1, 5, 
                                                                   settings['error_bar_thickness'],
                                                                   key=f"errthick_{gene}")
                        settings['error_bar_width'] = st.slider("Error Bar Width", 2, 10, 
                                                               settings['error_bar_width'],
                                                               key=f"errwidth_{gene}")
                    
                    settings['show_significance'] = st.checkbox("Significance Stars", 
                                                               settings['show_significance'],
                                                               key=f"sig_{gene}")
                    if settings['show_significance']:
                        settings['sig_position'] = st.selectbox("Star Position", 
                            ['outside', 'inside', 'auto'], key=f"sigpos_{gene}")
                
                with col3:
                    settings['show_grid'] = st.checkbox("Grid", settings['show_grid'], 
                                                       key=f"grid_{gene}")
                    settings['show_reference_line'] = st.checkbox("Reference Line", 
                                                                 settings['show_reference_line'],
                                                                 key=f"refline_{gene}")
                    if settings['show_reference_line']:
                        settings['reference_value'] = st.number_input("Reference Value", 
                            value=settings['reference_value'], key=f"refval_{gene}")
                    
                    settings['bar_border_width'] = st.slider("Bar Border Width", 0, 5, 
                                                            settings['bar_border_width'],
                                                            key=f"border_{gene}")
            
            # TAB: Colors
            with tab_colors:
                st.markdown("**Customize color for each condition**")
                
                # Initialize colors if needed
                if not settings['condition_colors']:
                    default_colors = px.colors.qualitative.Plotly
                    for i, cond in enumerate(settings['included_conditions']):
                        settings['condition_colors'][cond] = default_colors[i % len(default_colors)]
                
                cols = st.columns(3)
                for i, condition in enumerate(settings['included_conditions']):
                    with cols[i % 3]:
                        if condition not in settings['condition_colors']:
                            settings['condition_colors'][condition] = '#636EFA'
                        
                        settings['condition_colors'][condition] = st.color_picker(
                            condition,
                            settings['condition_colors'][condition],
                            key=f"color_{gene}_{condition}"
                        )
                
                col1, col2 = st.columns(2)
                with col1:
                    settings['sig_color'] = st.color_picker("Significance Color", 
                                                           settings['sig_color'],
                                                           key=f"sigcol_{gene}")
                with col2:
                    settings['grid_color'] = st.color_picker("Grid Color", 
                                                            settings['grid_color'],
                                                            key=f"gridcol_{gene}")
            
            # TAB: Labels
            with tab_labels:
                col1, col2 = st.columns(2)
                
                with col1:
                    settings['title'] = st.text_input("Title", settings['title'], 
                                                     key=f"title_{gene}")
                    settings['xlabel'] = st.text_input("X-axis Label", settings['xlabel'], 
                                                      key=f"xlabel_{gene}")
                    settings['ylabel'] = st.text_input("Y-axis Label", settings['ylabel'], 
                                                      key=f"ylabel_{gene}")
                
                with col2:
                    settings['title_size'] = st.slider("Title Size", 12, 36, settings['title_size'],
                                                      key=f"titsize_{gene}")
                    settings['axis_label_size'] = st.slider("Axis Label Size", 10, 24, 
                                                           settings['axis_label_size'],
                                                           key=f"axsize_{gene}")
                    settings['tick_font_size'] = st.slider("Tick Font Size", 8, 20, 
                                                          settings['tick_font_size'],
                                                          key=f"ticksize_{gene}")
                    settings['sig_font_size'] = st.slider("Significance Size", 10, 24, 
                                                         settings['sig_font_size'],
                                                         key=f"sigsize_{gene}")
                
                settings['x_tick_angle'] = st.slider("X-axis Tick Angle", -90, 90, 
                                                     settings['x_tick_angle'],
                                                     key=f"tickangle_{gene}")
            
            # TAB: Advanced
            with tab_advanced:
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    st.markdown("**Figure Dimensions**")
                    settings['figure_width'] = st.number_input("Width (px)", 400, 2000, 
                                                              settings['figure_width'],
                                                              key=f"width_{gene}")
                    settings['figure_height'] = st.number_input("Height (px)", 300, 1500, 
                                                               settings['figure_height'],
                                                               key=f"height_{gene}")
                
                with col2:
                    st.markdown("**Margins**")
                    settings['margin_left'] = st.number_input("Left", 0, 200, settings['margin_left'],
                                                             key=f"ml_{gene}")
                    settings['margin_right'] = st.number_input("Right", 0, 200, settings['margin_right'],
                                                              key=f"mr_{gene}")
                    settings['margin_top'] = st.number_input("Top", 0, 200, settings['margin_top'],
                                                            key=f"mt_{gene}")
                    settings['margin_bottom'] = st.number_input("Bottom", 0, 200, settings['margin_bottom'],
                                                               key=f"mb_{gene}")
                
                with col3:
                    st.markdown("**Y-axis Range**")
                    if st.checkbox(f"Custom Y range for {gene}", key=f"yrange_{gene}"):
                        y_min = st.number_input("Y min", value=0.0, key=f"ymin_{gene}")
                        y_max = st.number_input("Y max", value=3.0, key=f"ymax_{gene}")
                        settings['y_range'] = [y_min, y_max]
                    else:
                        settings['y_range'] = None
                    
                    st.markdown("**Background Colors**")
                    settings['plot_bgcolor'] = st.color_picker("Plot BG", settings['plot_bgcolor'],
                                                               key=f"plotbg_{gene}")
                    settings['paper_bgcolor'] = st.color_picker("Paper BG", settings['paper_bgcolor'],
                                                                key=f"paperbg_{gene}")
            
            # Generate graph
            st.markdown("### 📊 Live Preview")
            
            fig = GraphGenerator.create_gene_graph(
                gene_data,
                gene,
                settings,
                settings['included_conditions'],
                efficacy_config
            )
            
            st.plotly_chart(fig, use_container_width=True, key=f"preview_{gene}")
            
            if 'graphs' not in st.session_state:
                st.session_state.graphs = {}
            st.session_state.graphs[gene] = fig
            
            st.markdown("---")
    else:
        st.warning("⚠️ Run analysis first")

# ==================== TAB 5: EXPORT ====================
with tab5:
    st.header("Step 5: Export Everything")
    
    if st.session_state.processed_data:
        analysis_params = {
            'Date': datetime.now().strftime("%Y-%m-%d %H:%M"),
            'Efficacy_Type': st.session_state.selected_efficacy,
            'Housekeeping_Gene': st.session_state.hk_gene,
            'Baseline_Reference': baseline_condition if 'baseline_condition' in locals() else 'N/A',
            'Pvalue_Reference': pvalue_condition if 'pvalue_condition' in locals() else 'N/A',
            'Excluded_Wells': len(st.session_state.excluded_wells),
            'Excluded_Samples': len(st.session_state.excluded_samples),
            'Genes_Analyzed': len(st.session_state.processed_data)
        }
        
        st.subheader("📦 Complete Package")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("### 📊 Excel Report")
            st.caption("All calculations, raw data, and statistics")
            
            excel_data = export_to_excel(
                st.session_state.data,
                st.session_state.processed_data,
                analysis_params,
                st.session_state.sample_mapping
            )
            
            st.download_button(
                "📥 Download Excel",
                excel_data,
                f"qPCR_{st.session_state.selected_efficacy}_{datetime.now().strftime('%Y%m%d_%H%M')}.xlsx",
                "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                type="primary"
            )
        
        with col2:
            st.markdown("### 📈 All Graphs (HTML)")
            st.caption("Interactive graphs - works in PowerPoint!")
            
            if st.session_state.graphs:
                html_parts = [
                    "<html><head><title>qPCR Analysis</title></head><body>",
                    f"<h1>{st.session_state.selected_efficacy} Analysis</h1>",
                    f"<p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>"
                ]
                
                for gene, fig in st.session_state.graphs.items():
                    html_parts.append(f"<h2>{gene}</h2>")
                    html_parts.append(fig.to_html(include_plotlyjs='cdn', div_id=f"g_{gene}"))
                    html_parts.append("<hr>")
                
                html_parts.append("</body></html>")
                
                st.download_button(
                    "📥 Download All Graphs",
                    "\n".join(html_parts),
                    f"graphs_{st.session_state.selected_efficacy}_{datetime.now().strftime('%Y%m%d')}.html",
                    "text/html",
                    type="primary"
                )
        
        st.markdown("---")
        st.subheader("📋 Individual Files")
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.markdown("**Gene CSV Files**")
            for gene, df in st.session_state.processed_data.items():
                csv_buf = io.StringIO()
                df.drop(columns=['Target_Ct_Values'], errors='ignore').to_csv(csv_buf, index=False)
                st.download_button(
                    f"📥 {gene}.csv",
                    csv_buf.getvalue(),
                    f"{gene}_{datetime.now().strftime('%Y%m%d')}.csv",
                    "text/csv",
                    key=f"csv_{gene}"
                )
        
        with col2:
            st.markdown("**Gene HTML Graphs**")
            for gene, fig in st.session_state.graphs.items():
                html_buf = io.StringIO()
                fig.write_html(html_buf)
                st.download_button(
                    f"📥 {gene}.html",
                    html_buf.getvalue(),
                    f"{gene}_{datetime.now().strftime('%Y%m%d')}.html",
                    "text/html",
                    key=f"html_{gene}"
                )
        
        with col3:
            st.markdown("**Config Files**")
            config = {
                'params': analysis_params,
                'mapping': st.session_state.sample_mapping,
                'graph_settings': st.session_state.graph_settings,
                'included_samples': st.session_state.included_samples_for_analysis
            }
            st.download_button(
                "📥 Full Config",
                json.dumps(config, indent=2),
                f"config_{datetime.now().strftime('%Y%m%d_%H%M')}.json",
                "application/json"
            )
        
        with st.expander("💡 Export Tips"):
            st.markdown("""
            ### Publication-Ready Workflow
            1. **Excel**: Complete data for reviewers
            2. **HTML Graphs**: Open in browser → Right-click → Save as PNG/SVG
            3. **Config JSON**: Reproducibility proof
            
            ### For Presentations
            - Drag HTML files directly into PowerPoint/Google Slides
            - Interactive and zoomable!
            
            ### Browser Image Export
            - **Chrome/Edge**: Right-click graph → "Save image as..."
            - **Firefox**: Right-click → "Take Screenshot"
            - **Safari**: Right-click → "Export as Image"
            """)
        
        st.success("✅ All exports ready!")
    else:
        st.warning("⚠️ Complete analysis first")

# Footer
st.markdown("---")
st.markdown("""
<div style='text-align: center; color: #666;'>
    <p>🧬 qPCR Analysis Suite Ultimate v3.0</p>
    <p>Baseline reference • Custom sample selection • Full graph control • Publication-ready exports</p>
</div>
""", unsafe_allow_html=True)
