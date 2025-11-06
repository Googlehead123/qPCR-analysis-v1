import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from scipy import stats
import io
import json
from datetime import datetime
from typing import Dict, List, Tuple

# ==================== PAGE CONFIG ====================
st.set_page_config(page_title="qPCR Analysis Suite Pro", layout="wide", initial_sidebar_state="expanded")

# ==================== SESSION STATE INIT ====================
for key in ['data', 'processed_data', 'sample_mapping', 'analysis_templates', 'graphs', 
            'excluded_wells', 'excluded_samples', 'selected_efficacy', 'hk_gene']:
    if key not in st.session_state:
        st.session_state[key] = {} if key in ['sample_mapping', 'analysis_templates', 'graphs'] else (set() if 'excluded' in key else None)

# ==================== EFFICACY DATABASE ====================
EFFICACY_CONFIG = {
    '탄력': {
        'genes': ['COL1A1', 'ELN', 'FBN-1', 'FBN1'],
        'cell': 'HS68 fibroblast',
        'controls': {
            'negative': 'Non-treated',
            'positive': 'TGFb',
            'compare_to': 'negative'
        },
        'description': 'Elasticity - Non-treated vs TGFb (positive) vs Treatments'
    },
    '항노화': {
        'genes': ['COL1A1', 'COL1', 'MMP-1', 'MMP1'],
        'cell': 'HS68 fibroblast',
        'controls': {
            'baseline': 'Non-treated (No UV)',
            'negative': 'UVB only',
            'positive': 'UVB+TGFb',
            'compare_to': 'negative'
        },
        'description': 'Anti-aging - COL1↑ (recovery), MMP1↓ (inhibition) after UVB damage',
        'expected_direction': {'COL1A1': 'up', 'COL1': 'up', 'MMP-1': 'down', 'MMP1': 'down'}
    },
    '보습': {
        'genes': ['AQP3', 'HAS3'],
        'cell': 'HaCaT keratinocyte',
        'controls': {
            'negative': 'Non-treated',
            'positive': 'Retinoic acid',
            'compare_to': 'negative'
        },
        'description': 'Hydration - Non-treated vs Retinoic acid (positive) vs Treatments'
    },
    '장벽': {
        'genes': ['FLG', 'CLDN', 'IVL'],
        'cell': 'HaCaT keratinocyte',
        'controls': {
            'negative': 'Non-treated',
            'positive': 'Retinoic acid',
            'compare_to': 'negative'
        },
        'description': 'Barrier function - Non-treated vs Retinoic acid (positive) vs Treatments'
    },
    '표피증식': {
        'genes': ['KI67', 'PCNA'],
        'cell': 'HaCaT keratinocyte',
        'controls': {
            'negative': 'Non-treated',
            'positive': 'TGFb or FBS',
            'compare_to': 'negative'
        },
        'description': 'Proliferation - Non-treated vs TGFb/FBS (positive) vs Treatments'
    },
    '멜라닌억제': {
        'genes': ['MITF', 'TYR', 'Melanin'],
        'cell': 'B16F10 melanocyte',
        'controls': {
            'baseline': 'Non-treated',
            'negative': 'α-MSH only',
            'positive': 'α-MSH+Arbutin',
            'compare_to': 'negative'
        },
        'description': 'Melanin inhibition - α-MSH induced vs α-MSH+Arbutin (positive) vs α-MSH+Treatments',
        'expected_direction': {'MITF': 'down', 'TYR': 'down', 'Melanin': 'down'}
    },
    '진정': {
        'genes': ['IL1B', 'IL-1β', 'IL6', 'TNFA', 'TNFα'],
        'cell': 'HaCaT keratinocyte',
        'controls': {
            'baseline': 'Non-treated',
            'negative': 'IL4+PolyIC (Inflammation)',
            'positive': 'Inflammation+Dexamethasone',
            'compare_to': 'negative'
        },
        'description': 'Anti-inflammation - Reduce IL1β/IL6/TNFα (all should decrease)',
        'expected_direction': {'IL1B': 'down', 'IL-1β': 'down', 'IL6': 'down', 'TNFA': 'down', 'TNFα': 'down'}
    },
    '지질억제': {
        'genes': ['SREBPA', 'SREBPa', 'SREBPC', 'SREBPc', 'PPARY', 'PPARy'],
        'cell': 'SZ95 sebocyte',
        'controls': {
            'baseline': 'Non-treated',
            'negative': 'IGF only',
            'positive': 'IGF+Reference inhibitor',
            'compare_to': 'negative'
        },
        'description': 'Sebum inhibition - IGF induced vs IGF+Treatments',
        'expected_direction': {'SREBPA': 'down', 'SREBPa': 'down', 'SREBPC': 'down', 'SREBPc': 'down', 'PPARY': 'down', 'PPARy': 'down'}
    },
    '냉감': {
        'genes': ['TRPM8', 'CIRBP'],
        'cell': 'HaCaT keratinocyte',
        'controls': {
            'negative': 'Non-treated',
            'positive': 'Menthol',
            'compare_to': 'negative'
        },
        'description': 'Cooling effect - Non-treated vs Menthol (positive) vs Treatments'
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
                return ('format2' if 'Cт' in row_str or 'ΔCт' in row_str else 'format1'), idx
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
        
        parsed = pd.DataFrame({
            'Well': df['Well'],
            'Sample': df['Sample Name'],
            'Target': df['Target Name'],
            'CT': pd.to_numeric(df['Cт'], errors='coerce')
        })
        
        return parsed.dropna(subset=['CT']).query('Sample.notna() & Target.notna()')
    
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
    def calculate_ddct(data: pd.DataFrame, hk_gene: str, ref_sample: str, compare_sample: str, 
                       excluded_wells: set, excluded_samples: set, sample_mapping: dict) -> pd.DataFrame:
        """Gene-by-gene ΔΔCt calculation with housekeeping normalization"""
        
        # Filter data
        data = data[~data['Well'].isin(excluded_wells) & ~data['Sample'].isin(excluded_samples)].copy()
        
        # Apply sample name mapping
        data['Condition'] = data['Sample'].map(lambda x: sample_mapping.get(x, {}).get('condition', x))
        data['Group'] = data['Sample'].map(lambda x: sample_mapping.get(x, {}).get('group', 'Treatment'))
        
        results = []
        
        # Process each target gene separately (exclude housekeeping)
        for target in data['Target'].unique():
            if target.upper() in [hk_gene.upper(), 'ACTIN', 'B-ACTIN', 'GAPDH', 'ACTB']:
                continue
            
            target_data = data[data['Target'] == target]
            
            for condition in target_data['Condition'].unique():
                cond_data = target_data[target_data['Condition'] == condition]
                
                # Get housekeeping Ct for this condition
                hk_data = data[(data['Condition'] == condition) & (data['Target'] == hk_gene)]
                if len(hk_data) == 0:
                    continue
                
                # Calculate ΔCt = Target_Ct - HK_Ct
                target_ct_values = cond_data['CT'].values
                hk_ct_values = hk_data['CT'].values
                
                # Use mean if multiple replicates
                target_ct_mean = target_ct_values.mean()
                hk_ct_mean = hk_ct_values.mean()
                delta_ct = target_ct_mean - hk_ct_mean
                
                # Get reference ΔCt (compare_sample)
                ref_target = target_data[target_data['Condition'] == compare_sample]
                ref_hk = data[(data['Condition'] == compare_sample) & (data['Target'] == hk_gene)]
                
                if len(ref_target) > 0 and len(ref_hk) > 0:
                    ref_delta_ct = ref_target['CT'].mean() - ref_hk['CT'].mean()
                else:
                    ref_delta_ct = 0
                
                # ΔΔCt and relative expression
                ddct = delta_ct - ref_delta_ct
                rel_expr = 2 ** (-ddct)
                
                # Error calculation
                ct_sd = target_ct_values.std() if len(target_ct_values) > 1 else 0
                sem = ct_sd / np.sqrt(len(target_ct_values)) if len(target_ct_values) > 1 else 0
                
                # Get original sample name and group
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
                    'HK_Ct_Mean': hk_ct_mean,
                    'Delta_Ct': delta_ct,
                    'Delta_Delta_Ct': ddct,
                    'Relative_Expression': rel_expr,
                    'SEM': sem,
                    'Fold_Change': rel_expr
                })
        
        return pd.DataFrame(results)
    
    @staticmethod
    def calculate_statistics(data: pd.DataFrame, compare_condition: str) -> pd.DataFrame:
        """Calculate p-values comparing each condition to comparison control"""
        results = data.copy()
        results['p_value'] = np.nan
        results['significance'] = ''
        
        for target in results['Target'].unique():
            target_data = results[results['Target'] == target]
            control_values = target_data[target_data['Condition'] == compare_condition]['Relative_Expression'].values
            
            for condition in target_data['Condition'].unique():
                if condition == compare_condition:
                    continue
                
                sample_values = target_data[target_data['Condition'] == condition]['Relative_Expression'].values
                
                if len(control_values) >= 1 and len(sample_values) >= 1:
                    # Use one-sample t-test if only one control value
                    if len(control_values) == 1 and len(sample_values) > 1:
                        t_stat, p_val = stats.ttest_1samp(sample_values, control_values[0])
                    elif len(sample_values) == 1 and len(control_values) > 1:
                        t_stat, p_val = stats.ttest_1samp(control_values, sample_values[0])
                    else:
                        t_stat, p_val = stats.ttest_ind(control_values, sample_values)
                    
                    idx = (results['Target'] == target) & (results['Condition'] == condition)
                    results.loc[idx, 'p_value'] = p_val
                    
                    if p_val < 0.001:
                        results.loc[idx, 'significance'] = '***'
                    elif p_val < 0.01:
                        results.loc[idx, 'significance'] = '**'
                    elif p_val < 0.05:
                        results.loc[idx, 'significance'] = '*'
        
        return results

# ==================== GRAPH GENERATOR ====================
class GraphGenerator:
    @staticmethod
    def create_gene_graph(data: pd.DataFrame, gene: str, settings: dict, efficacy_config: dict = None) -> go.Figure:
        """Create individual graph for each gene"""
        gene_data = data[data['Target'] == gene].copy()
        
        # Sort by group order: baseline → negative → positive → treatments
        group_order = ['Baseline', 'Negative Control', 'Positive Control', 'Treatment']
        gene_data['group_sort'] = gene_data['Group'].apply(
            lambda x: group_order.index(x) if x in group_order else 999
        )
        gene_data = gene_data.sort_values(['group_sort', 'Condition'])
        
        fig = go.Figure()
        
        # Color mapping by group
        color_map = {
            'Baseline': '#808080',
            'Negative Control': '#FF6B6B',
            'Positive Control': '#4ECDC4',
            'Treatment': settings.get('bar_colors', {}).get(gene, '#95E1D3')
        }
        
        fig.add_trace(go.Bar(
            x=gene_data['Condition'],
            y=gene_data['Relative_Expression'],
            error_y=dict(
                type='data',
                array=gene_data['SEM'] * 1.96,
                visible=settings.get('show_error', True)
            ),
            text=gene_data['significance'] if settings.get('show_significance', True) else None,
            textposition='outside',
            textfont=dict(size=settings.get('sig_font_size', 16)),
            marker=dict(color=[color_map.get(g, '#95E1D3') for g in gene_data['Group']]),
            showlegend=False
        ))
        
        # Add expected direction indicator
        if efficacy_config and 'expected_direction' in efficacy_config:
            direction = efficacy_config['expected_direction'].get(gene, '')
            direction_text = '↑ Expected increase' if direction == 'up' else '↓ Expected decrease' if direction == 'down' else ''
            if direction_text:
                fig.add_annotation(
                    text=direction_text,
                    xref='paper', yref='paper',
                    x=0.02, y=0.98,
                    showarrow=False,
                    font=dict(size=12, color='red'),
                    align='left'
                )
        
        fig.update_layout(
            title=dict(text=f"{gene} Expression", font=dict(size=settings.get('title_size', 20))),
            xaxis_title=settings.get('xlabel', 'Condition'),
            yaxis_title=settings.get('ylabel', 'Fold Change (Relative to Control)'),
            template=settings.get('color_scheme', 'plotly_white'),
            font=dict(size=settings.get('font_size', 14)),
            height=settings.get('figure_height', 600),
            width=settings.get('figure_width', 1000),
            xaxis=dict(showgrid=settings.get('show_grid', True)),
            yaxis=dict(showgrid=settings.get('show_grid', True))
        )
        
        return fig

# ==================== EXPORT FUNCTIONS ====================
def export_to_excel(raw_data: pd.DataFrame, processed_data: Dict[str, pd.DataFrame], 
                   params: dict, mapping: dict) -> bytes:
    """Export comprehensive Excel with gene-by-gene sheets"""
    output = io.BytesIO()
    
    with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
        # Parameters sheet
        pd.DataFrame([params]).to_excel(writer, sheet_name='Analysis_Parameters', index=False)
        
        # Sample mapping sheet
        pd.DataFrame([{'Original': k, **v} for k, v in mapping.items()]).to_excel(
            writer, sheet_name='Sample_Mapping', index=False
        )
        
        # Raw data
        raw_data.to_excel(writer, sheet_name='Raw_Data', index=False)
        
        # Gene-by-gene calculations
        for gene, gene_data in processed_data.items():
            sheet_name = f"{gene}_Analysis"[:31]  # Excel sheet name limit
            gene_data.to_excel(writer, sheet_name=sheet_name, index=False)
        
        # Summary sheet
        if processed_data:
            all_data = pd.concat(processed_data.values(), ignore_index=True)
            summary = all_data.groupby(['Target', 'Group']).agg({
                'Relative_Expression': ['mean', 'std', 'count'],
                'p_value': 'min'
            }).round(4)
            summary.to_excel(writer, sheet_name='Summary')
    
    return output.getvalue()

# ==================== UI ====================
st.title("🧬 qPCR Analysis Suite Pro")
st.markdown("**Gene-by-gene analysis with efficacy-specific workflows**")

# Sidebar
with st.sidebar:
    st.header("💬 Quick Guide")
    st.markdown("""
    1. **Upload** CSV files
    2. **Select** efficacy type
    3. **Map** samples to conditions
    4. **Exclude** unwanted samples/wells
    5. **Analyze** gene-by-gene
    6. **Customize** individual graphs
    7. **Export** all results
    """)
    
    st.subheader("📋 Templates")
    template_name = st.text_input("Save analysis as:")
    if st.button("💾 Save") and st.session_state.sample_mapping:
        st.session_state.analysis_templates[template_name] = {
            'mapping': st.session_state.sample_mapping.copy(),
            'efficacy': st.session_state.selected_efficacy,
            'timestamp': datetime.now().isoformat()
        }
        st.success(f"✅ Saved '{template_name}'")
    
    if st.session_state.analysis_templates:
        load = st.selectbox("Load:", [""] + list(st.session_state.analysis_templates.keys()))
        if load:
            template = st.session_state.analysis_templates[load]
            st.session_state.sample_mapping = template['mapping']
            st.session_state.selected_efficacy = template.get('efficacy')
            st.info(f"✅ Loaded '{load}'")

# Main tabs
tab1, tab2, tab3, tab4, tab5 = st.tabs(["📁 Upload & Filter", "🗺️ Sample Mapping", "🔬 Analysis", "📊 Graphs", "📤 Export"])

# ==================== TAB 1: UPLOAD & FILTER ====================
with tab1:
    st.header("Step 1: Upload & Filter Data")
    
    uploaded_files = st.file_uploader("Upload qPCR CSV files", type=['csv'], accept_multiple_files=True)
    
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
            
            # Detect housekeeping gene
            hk_genes = [g for g in st.session_state.data['Target'].unique() 
                       if g.upper() in ['ACTIN', 'B-ACTIN', 'GAPDH', 'ACTB']]
            if hk_genes:
                st.session_state.hk_gene = st.selectbox("🔬 Housekeeping Gene", hk_genes, key='hk_select')
                col4.metric("HK Gene", st.session_state.hk_gene)
            
            # Sample filter
            st.subheader("🎯 Filter Samples")
            st.markdown("*Select samples to INCLUDE in analysis*")
            
            all_samples = sorted(st.session_state.data['Sample'].unique())
            selected_samples = st.multiselect(
                "Include these samples:",
                all_samples,
                default=all_samples,
                help="Deselect samples you want to exclude"
            )
            st.session_state.excluded_samples = set(all_samples) - set(selected_samples)
            
            if st.session_state.excluded_samples:
                st.warning(f"⚠️ Excluding {len(st.session_state.excluded_samples)} samples: {', '.join(list(st.session_state.excluded_samples)[:5])}")
            
            # Gene filter
            st.subheader("🧬 Filter Genes")
            all_genes = sorted(st.session_state.data['Target'].unique())
            selected_genes = st.multiselect(
                "Include these genes:",
                all_genes,
                default=all_genes,
                help="Deselect genes you want to exclude from analysis"
            )
            
            # Filter data display
            display_data = st.session_state.data[
                st.session_state.data['Sample'].isin(selected_samples) &
                st.session_state.data['Target'].isin(selected_genes)
            ]
            
            st.subheader("📊 Filtered Data Preview")
            st.dataframe(display_data.head(50), height=300)
            
            # Outlier detection
            st.subheader("🚩 Outlier Detection")
            for target in selected_genes:
                target_data = display_data[display_data['Target'] == target]
                q1, q3 = target_data['CT'].quantile([0.25, 0.75])
                iqr = q3 - q1
                outliers = target_data[(target_data['CT'] < q1 - 1.5*iqr) | (target_data['CT'] > q3 + 1.5*iqr)]
                
                if len(outliers) > 0:
                    with st.expander(f"⚠️ {target}: {len(outliers)} outliers"):
                        st.dataframe(outliers[['Well', 'Sample', 'CT']])
                        if st.checkbox(f"Exclude these wells", key=f"out_{target}"):
                            st.session_state.excluded_wells.update(outliers['Well'].tolist())

# ==================== TAB 2: SAMPLE MAPPING ====================
with tab2:
    st.header("Step 2: Map Samples to Conditions")
    
    if st.session_state.data is not None:
        # Efficacy type selection
        detected_genes = set(st.session_state.data['Target'].unique())
        suggested = None
        for eff, cfg in EFFICACY_CONFIG.items():
            if any(g in detected_genes for g in cfg['genes']):
                suggested = eff
                break
        
        efficacy = st.selectbox(
            "🎯 Efficacy Test Type",
            list(EFFICACY_CONFIG.keys()),
            index=list(EFFICACY_CONFIG.keys()).index(suggested) if suggested else 0
        )
        st.session_state.selected_efficacy = efficacy
        
        config = EFFICACY_CONFIG[efficacy]
        st.info(f"**{config['description']}**")
        st.caption(f"Cell line: {config['cell']} | Genes: {', '.join(config['genes'])}")
        
        # Show control structure
        with st.expander("📋 Control Structure for this Test"):
            for ctrl_type, ctrl_name in config['controls'].items():
                st.markdown(f"- **{ctrl_type.title()}**: {ctrl_name}")
        
        # Sample mapping interface
        st.subheader("🗺️ Sample Mapping")
        
        samples = [s for s in sorted(st.session_state.data['Sample'].unique()) 
                  if s not in st.session_state.excluded_samples]
        
        # Group type options based on efficacy
        group_types = ['Negative Control', 'Positive Control', 'Treatment']
        if 'baseline' in config['controls']:
            group_types.insert(0, 'Baseline')
        
        for sample in samples:
            if sample not in st.session_state.sample_mapping:
                st.session_state.sample_mapping[sample] = {
                    'condition': sample,
                    'group': 'Treatment',
                    'concentration': ''
                }
            
            col1, col2, col3, col4 = st.columns([1, 3, 2, 2])
            
            with col1:
                st.text(sample)
            with col2:
                condition = st.text_input(
                    "Condition name",
                    value=st.session_state.sample_mapping[sample]['condition'],
                    key=f"cond_{sample}",
                    label_visibility="collapsed"
                )
                st.session_state.sample_mapping[sample]['condition'] = condition
            with col3:
                group = st.selectbox(
                    "Group",
                    group_types,
                    index=group_types.index(st.session_state.sample_mapping[sample]['group']) 
                          if st.session_state.sample_mapping[sample]['group'] in group_types else 0,
                    key=f"grp_{sample}",
                    label_visibility="collapsed"
                )
                st.session_state.sample_mapping[sample]['group'] = group
            with col4:
                conc = st.text_input(
                    "Concentration",
                    value=st.session_state.sample_mapping[sample]['concentration'],
                    key=f"conc_{sample}",
                    label_visibility="collapsed"
                )
                st.session_state.sample_mapping[sample]['concentration'] = conc
        
        # Summary
        st.subheader("📊 Mapping Summary")
        mapping_df = pd.DataFrame([{'Original': k, **v} for k, v in st.session_state.sample_mapping.items()])
        st.dataframe(mapping_df, use_container_width=True)
    else:
        st.warning("⚠️ Upload data first")

# ==================== TAB 3: ANALYSIS ====================
with tab3:
    st.header("Step 3: Gene-by-Gene Analysis")
    
    if st.session_state.data is not None and st.session_state.sample_mapping and st.session_state.hk_gene:
        config = EFFICACY_CONFIG.get(st.session_state.selected_efficacy, {})
        
        col1, col2 = st.columns(2)
        
        with col1:
            # Select reference sample (what treatments are relative to)
            ref_samples = [k for k, v in st.session_state.sample_mapping.items() 
                          if v['group'] == 'Negative Control']
            if not ref_samples and 'baseline' in config.get('controls', {}):
                ref_samples = [k for k, v in st.session_state.sample_mapping.items() 
                              if v['group'] == 'Baseline']
            
            ref_sample = st.selectbox(
                "Reference Sample (Fold change = 1.0)",
                ref_samples if ref_samples else list(st.session_state.sample_mapping.keys())
            )
            ref_condition = st.session_state.sample_mapping[ref_sample]['condition']
        
        with col2:
            # Select comparison control (for p-value calculation)
            compare_samples = [k for k, v in st.session_state.sample_mapping.items() 
                             if v['group'] == 'Negative Control']
            compare_sample = st.selectbox(
                "Compare To (for p-values)",
                compare_samples if compare_samples else list(st.session_state.sample_mapping.keys())
            )
            compare_condition = st.session_state.sample_mapping[compare_sample]['condition']
        
        st.info(f"📊 Analysis setup: All fold changes relative to **{ref_condition}**, p-values vs **{compare_condition}**")
        
        if st.button("🔬 Run Gene-by-Gene Analysis", type="primary"):
            with st.spinner("Analyzing each gene separately..."):
                # Calculate ΔΔCt
                processed = AnalysisEngine.calculate_ddct(
                    st.session_state.data,
                    st.session_state.hk_gene,
                    ref_condition,
                    compare_condition,
                    st.session_state.excluded_wells,
                    st.session_state.excluded_samples,
                    st.session_state.sample_mapping
                )
                
                if processed is not None and len(processed) > 0:
                    # Calculate statistics
                    processed = AnalysisEngine.calculate_statistics(processed, compare_condition)
                    
                    # Split by gene
                    gene_data = {}
                    for gene in processed['Target'].unique():
                        gene_df = processed[processed['Target'] == gene].copy()
                        gene_data[gene] = gene_df
                    
                    st.session_state.processed_data = gene_data
                    st.success(f"✅ Analysis complete! {len(gene_data)} genes analyzed")
        
        # Display results
        if st.session_state.processed_data:
            st.subheader("📊 Analysis Results by Gene")
            
            # Overall summary
            all_results = pd.concat(st.session_state.processed_data.values(), ignore_index=True)
            col1, col2, col3, col4 = st.columns(4)
            col1.metric("Genes Analyzed", len(st.session_state.processed_data))
            col2.metric("Conditions", all_results['Condition'].nunique())
            sig_count = (all_results['p_value'] < 0.05).sum()
            col3.metric("Significant (p<0.05)", sig_count)
            col4.metric("Total Comparisons", len(all_results))
            
            # Gene-by-gene tables
            for gene, gene_df in st.session_state.processed_data.items():
                with st.expander(f"🧬 {gene} - Detailed Results", expanded=True):
                    # Show expected direction if available
                    config = EFFICACY_CONFIG.get(st.session_state.selected_efficacy, {})
                    if 'expected_direction' in config and gene in config['expected_direction']:
                        direction = config['expected_direction'][gene]
                        st.caption(f"Expected: {'↑ Increase' if direction == 'up' else '↓ Decrease'}")
                    
                    display_cols = ['Condition', 'Group', 'Fold_Change', 'p_value', 'significance', 
                                  'n_replicates', 'Target_Ct_Mean', 'SEM']
                    
                    styled_df = gene_df[display_cols].style.background_gradient(
                        subset=['Fold_Change'], cmap='RdYlGn', vmin=0, vmax=3
                    ).format({
                        'Fold_Change': '{:.3f}',
                        'p_value': '{:.4f}',
                        'Target_Ct_Mean': '{:.2f}',
                        'SEM': '{:.3f}'
                    })
                    
                    st.dataframe(styled_df, use_container_width=True)
    else:
        st.warning("⚠️ Complete previous steps first")

# ==================== TAB 4: GRAPHS ====================
with tab4:
    st.header("Step 4: Individual Gene Graphs")
    
    if st.session_state.processed_data:
        # Initialize graph settings
        if 'graph_settings' not in st.session_state:
            st.session_state.graph_settings = {
                'title_size': 20, 'font_size': 14, 'sig_font_size': 16,
                'figure_width': 1000, 'figure_height': 600,
                'color_scheme': 'plotly_white', 'show_error': True,
                'show_significance': True, 'show_grid': True,
                'xlabel': 'Condition', 'ylabel': 'Fold Change (Relative to Control)',
                'bar_colors': {}
            }
        
        # Global settings
        st.subheader("⚙️ Global Graph Settings")
        
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.session_state.graph_settings['font_size'] = st.slider("Font Size", 10, 24, 14)
            st.session_state.graph_settings['title_size'] = st.slider("Title Size", 14, 32, 20)
        with col2:
            st.session_state.graph_settings['figure_width'] = st.slider("Width (px)", 600, 1600, 1000)
            st.session_state.graph_settings['figure_height'] = st.slider("Height (px)", 400, 1200, 600)
        with col3:
            st.session_state.graph_settings['color_scheme'] = st.selectbox(
                "Theme", ['plotly_white', 'plotly', 'plotly_dark', 'seaborn', 'simple_white', 'presentation']
            )
        with col4:
            st.session_state.graph_settings['show_error'] = st.checkbox("Error Bars", True)
            st.session_state.graph_settings['show_significance'] = st.checkbox("Significance", True)
            st.session_state.graph_settings['show_grid'] = st.checkbox("Grid", True)
        
        # Quick presets
        col1, col2, col3, col4 = st.columns(4)
        if col1.button("📄 Publication"):
            st.session_state.graph_settings.update({
                'color_scheme': 'simple_white', 'font_size': 14, 'title_size': 18,
                'figure_width': 1000, 'figure_height': 600
            })
            st.rerun()
        if col2.button("🎨 Presentation"):
            st.session_state.graph_settings.update({
                'color_scheme': 'presentation', 'font_size': 18, 'title_size': 24,
                'figure_width': 1200, 'figure_height': 700
            })
            st.rerun()
        if col3.button("🌙 Dark Mode"):
            st.session_state.graph_settings.update({'color_scheme': 'plotly_dark'})
            st.rerun()
        if col4.button("🔄 Reset"):
            del st.session_state.graph_settings
            st.rerun()
        
        st.markdown("---")
        
        # Generate graphs for each gene
        st.subheader("📊 Gene-Specific Graphs")
        
        efficacy_config = EFFICACY_CONFIG.get(st.session_state.selected_efficacy, {})
        
        # Create graphs
        if 'graphs' not in st.session_state:
            st.session_state.graphs = {}
        
        for gene in st.session_state.processed_data.keys():
            st.markdown(f"### 🧬 {gene}")
            
            # Gene-specific color picker
            col1, col2 = st.columns([1, 4])
            with col1:
                if gene not in st.session_state.graph_settings['bar_colors']:
                    default_colors = px.colors.qualitative.Plotly
                    idx = list(st.session_state.processed_data.keys()).index(gene)
                    st.session_state.graph_settings['bar_colors'][gene] = default_colors[idx % len(default_colors)]
                
                st.session_state.graph_settings['bar_colors'][gene] = st.color_picker(
                    f"{gene} color",
                    st.session_state.graph_settings['bar_colors'][gene],
                    key=f"color_{gene}"
                )
            
            # Generate graph
            fig = GraphGenerator.create_gene_graph(
                st.session_state.processed_data[gene],
                gene,
                st.session_state.graph_settings,
                efficacy_config
            )
            
            st.plotly_chart(fig, use_container_width=True, key=f"graph_{gene}")
            st.session_state.graphs[gene] = fig
            
            st.markdown("---")
    else:
        st.warning("⚠️ Run analysis first")

# ==================== TAB 5: EXPORT ====================
with tab5:
    st.header("Step 5: Export All Results")
    
    if st.session_state.processed_data:
        st.subheader("📦 Download Options")
        
        # Prepare export parameters
        analysis_params = {
            'Date': datetime.now().strftime("%Y-%m-%d %H:%M"),
            'Efficacy_Type': st.session_state.selected_efficacy,
            'Housekeeping_Gene': st.session_state.hk_gene,
            'Reference_Sample': ref_condition if 'ref_condition' in locals() else 'N/A',
            'Compare_To': compare_condition if 'compare_condition' in locals() else 'N/A',
            'Excluded_Wells': len(st.session_state.excluded_wells),
            'Excluded_Samples': len(st.session_state.excluded_samples),
            'Genes_Analyzed': len(st.session_state.processed_data)
        }
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("### 📊 Complete Excel Report")
            st.caption("Includes: Parameters, Mapping, Raw Data, Gene-by-Gene Calculations, Summary")
            
            excel_data = export_to_excel(
                st.session_state.data,
                st.session_state.processed_data,
                analysis_params,
                st.session_state.sample_mapping
            )
            
            st.download_button(
                label="📥 Download Excel Report",
                data=excel_data,
                file_name=f"qPCR_{st.session_state.selected_efficacy}_{datetime.now().strftime('%Y%m%d_%H%M')}.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                type="primary"
            )
        
        with col2:
            st.markdown("### 📈 All Graphs (HTML)")
            st.caption("Interactive graphs for all genes in one file")
            
            if st.session_state.graphs:
                # Create combined HTML
                html_parts = ["<html><head><title>qPCR Analysis Graphs</title></head><body>"]
                html_parts.append(f"<h1>{st.session_state.selected_efficacy} Analysis</h1>")
                html_parts.append(f"<p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>")
                
                for gene, fig in st.session_state.graphs.items():
                    html_parts.append(f"<h2>{gene}</h2>")
                    html_parts.append(fig.to_html(include_plotlyjs='cdn', div_id=f"graph_{gene}"))
                    html_parts.append("<hr>")
                
                html_parts.append("</body></html>")
                combined_html = "\n".join(html_parts)
                
                st.download_button(
                    label="📥 Download All Graphs (HTML)",
                    data=combined_html,
                    file_name=f"qPCR_graphs_{st.session_state.selected_efficacy}_{datetime.now().strftime('%Y%m%d_%H%M')}.html",
                    mime="text/html",
                    type="primary"
                )
        
        st.markdown("---")
        st.subheader("📋 Individual Downloads")
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            # CSV export per gene
            st.markdown("**Gene-by-Gene CSV**")
            for gene, gene_df in st.session_state.processed_data.items():
                csv_buffer = io.StringIO()
                gene_df.to_csv(csv_buffer, index=False)
                
                st.download_button(
                    label=f"📥 {gene}.csv",
                    data=csv_buffer.getvalue(),
                    file_name=f"{gene}_calculations_{datetime.now().strftime('%Y%m%d')}.csv",
                    mime="text/csv",
                    key=f"csv_{gene}"
                )
        
        with col2:
            # Individual graph HTML
            st.markdown("**Individual Graph HTML**")
            for gene, fig in st.session_state.graphs.items():
                html_buffer = io.StringIO()
                fig.write_html(html_buffer)
                
                st.download_button(
                    label=f"📥 {gene}.html",
                    data=html_buffer.getvalue(),
                    file_name=f"{gene}_graph_{datetime.now().strftime('%Y%m%d')}.html",
                    mime="text/html",
                    key=f"html_{gene}"
                )
        
        with col3:
            # Configuration export
            st.markdown("**Reproducibility Files**")
            
            # Analysis config
            config_data = {
                'analysis_params': analysis_params,
                'sample_mapping': st.session_state.sample_mapping,
                'graph_settings': st.session_state.graph_settings,
                'excluded_wells': list(st.session_state.excluded_wells),
                'excluded_samples': list(st.session_state.excluded_samples)
            }
            
            st.download_button(
                label="📥 Analysis Config",
                data=json.dumps(config_data, indent=2),
                file_name=f"config_{datetime.now().strftime('%Y%m%d_%H%M')}.json",
                mime="application/json"
            )
            
            # Graph preset
            st.download_button(
                label="📥 Graph Preset",
                data=json.dumps(st.session_state.graph_settings, indent=2),
                file_name=f"graph_preset_{datetime.now().strftime('%Y%m%d')}.json",
                mime="application/json"
            )
        
        # Export tips
        with st.expander("💡 Export Guide"):
            st.markdown(f"""
            ### For {st.session_state.selected_efficacy} Analysis
            
            **Complete Package (Recommended)**
            - ✅ Excel Report: All calculations, statistics, and raw data
            - ✅ All Graphs HTML: Interactive figures for all {len(st.session_state.graphs)} genes
            - ✅ Analysis Config: For reproducibility and audit trail
            
            **For Publications**
            1. Download Excel → Reviewer can verify calculations
            2. Download individual HTML → Open in browser → Right-click → Save as PNG/SVG
            3. Or screenshot from browser for manuscripts
            
            **For Presentations**
            - Drag HTML files directly into PowerPoint/Google Slides
            - Interactive graphs work in presentations!
            
            **For Patents/IP**
            - Excel: Complete audit trail with timestamps
            - Config JSON: Reproducibility proof
            - All Graphs HTML: Visual evidence
            
            **Gene-by-Gene Files**
            - Useful for sharing specific gene results
            - CSV for data analysis in other tools
            - HTML for interactive sharing
            """)
        
        st.success("✅ All export options ready!")
    else:
        st.warning("⚠️ Complete analysis first")

# ==================== FOOTER ====================
st.markdown("---")
st.markdown("""
<div style='text-align: center; color: #666;'>
    <p>🧬 qPCR Analysis Suite Pro v2.0 | Gene-by-gene analysis with efficacy-specific workflows</p>
    <p>Housekeeping normalization • Individual gene graphs • Comprehensive statistics</p>
</div>
""", unsafe_allow_html=True)
