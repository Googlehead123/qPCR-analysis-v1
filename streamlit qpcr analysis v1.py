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
                
                # Get reference ΔCt (ref_sample)
                ref_target = target_data[target_data['Condition'] == ref_sample]
                ref_hk = data[(data['Condition'] == ref_sample) & (data['Target'] == hk_gene)]

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
    def calculate_statistics(processed: pd.DataFrame, compare_condition: str,
                            raw_data: pd.DataFrame = None,
                            hk_gene: str = None,
                            sample_mapping: dict = None) -> pd.DataFrame:
        """Two-tailed Welch's t-test comparing each condition to compare_condition"""
        
        # Use session_state fallbacks
        raw_data = raw_data if raw_data is not None else st.session_state.get("data")
        hk_gene = hk_gene if hk_gene is not None else st.session_state.get("hk_gene")
        sample_mapping = sample_mapping if sample_mapping is not None else st.session_state.get("sample_mapping", {})
        
        if raw_data is None or hk_gene is None:
            return processed
        
        results = processed.copy()
        results["p_value"] = np.nan
        results["significance"] = ""
        
        for target in results["Target"].unique():
            if pd.isna(target):
                continue
            
            # Map conditions
            t_rows = raw_data[raw_data["Target"] == target].copy()
            if t_rows.empty:
                continue
            
            t_rows["Condition"] = t_rows["Sample"].map(
                lambda s: sample_mapping.get(s, {}).get("condition", s)
            )
            
            hk_rows = raw_data[raw_data["Target"] == hk_gene].copy()
            hk_rows["Condition"] = hk_rows["Sample"].map(
                lambda s: sample_mapping.get(s, {}).get("condition", s)
            )
            hk_means = hk_rows.groupby("Condition")["CT"].mean().to_dict()
            
            # Calculate relative expression per condition
            rel_expr = {}
            for cond, grp in t_rows.groupby("Condition"):
                hk_mean = hk_means.get(cond, np.nan)
                if np.isnan(hk_mean):
                    continue
                rel_expr[cond] = 2 ** (-(grp["CT"].values - hk_mean))
            
            # Reference condition
            ref_vals = rel_expr.get(compare_condition, np.array([]))
            if ref_vals.size < 1:
                continue
            
            # Compare each condition to reference
            for cond, vals in rel_expr.items():
                if cond == compare_condition or vals.size == 0:
                    continue
                
                try:
                    # Choose appropriate test
                    if ref_vals.size >= 2 and vals.size >= 2:
                        _, p_val = stats.ttest_ind(ref_vals, vals, equal_var=False)
                    elif vals.size == 1 and ref_vals.size >= 2:
                        _, p_val = stats.ttest_1samp(ref_vals, vals[0])
                    elif ref_vals.size == 1 and vals.size >= 2:
                        _, p_val = stats.ttest_1samp(vals, ref_vals[0])
                    else:
                        p_val = np.nan
                except:
                    p_val = np.nan
                
                # Annotate results
                mask = (results["Target"] == target) & (results["Condition"] == cond)
                results.loc[mask, "p_value"] = p_val
                
                if not np.isnan(p_val):
                    if p_val < 0.001:
                        results.loc[mask, "significance"] = "***"
                    elif p_val < 0.01:
                        results.loc[mask, "significance"] = "**"
                    elif p_val < 0.05:
                        results.loc[mask, "significance"] = "*"
        
        return results
    # ------------------------------------------------------------
    # Helper function: Run full qPCR ΔΔCt + statistics pipeline
    # ------------------------------------------------------------
    def run_full_analysis(ref_sample_key: str, compare_sample_key: str):
        """
        Run ΔΔCt + statistical analysis and store results in st.session_state.
        Produces st.session_state.processed_data = {gene: DataFrame}.
        """
        try:
            data = st.session_state.get("data")
            mapping = st.session_state.get("sample_mapping", {})
            hk_gene = st.session_state.get("hk_gene")

            if data is None:
                st.error("❌ No raw data loaded.")
                return False
            if not mapping:
                st.error("❌ Sample mapping not found.")
                return False
            if not hk_gene:
                st.error("❌ Housekeeping gene not selected.")
                return False

            ref_condition = mapping.get(ref_sample_key, {}).get("condition", ref_sample_key)
            cmp_condition = mapping.get(compare_sample_key, {}).get("condition", compare_sample_key)

            with st.spinner(f"Running full analysis using reference '{ref_condition}' and comparison '{cmp_condition}'..."):
                # --- ΔΔCt calculation ---
                processed_df = AnalysisEngine.calculate_ddct(
                    data,
                    hk_gene,
                    ref_condition,
                    cmp_condition,
                    st.session_state.get("excluded_wells", set()),
                    st.session_state.get("excluded_samples", set()),
                    mapping,
                )

                if processed_df is None or processed_df.empty:
                    st.warning("⚠️ No ΔΔCt results produced. Check mapping and housekeeping gene.")
                    return False

                # --- Statistical test ---
                try:
                    processed_with_stats = AnalysisEngine.calculate_statistics(
                        processed_df,
                        cmp_condition,
                        raw_data=data,
                        hk_gene=hk_gene,
                        sample_mapping=mapping,
                    )
                except TypeError:
                    # fallback for simpler signature
                    processed_with_stats = AnalysisEngine.calculate_statistics(processed_df, cmp_condition)

                # --- Organize data for graphs ---
                gene_dict = {}
                if "Target" in processed_with_stats.columns:
                    for gene in processed_with_stats["Target"].unique():
                        gene_df = processed_with_stats[processed_with_stats["Target"] == gene].copy()
                        gene_dict[gene] = gene_df.reset_index(drop=True)
                else:
                    gene_dict = {"results": processed_with_stats.reset_index(drop=True)}

                st.session_state.processed_data = gene_dict

            st.success("✅ Full analysis complete. Go to the Graphs tab to visualize results.")
            return True

        except Exception as e:
            st.error(f"❌ Analysis failed: {e}")
            return False

# ==================== GRAPH GENERATOR ====================
class GraphGenerator:
    @staticmethod
    def create_gene_graph(
        data: pd.DataFrame,
        gene: str,
        settings: dict,
        efficacy_config: dict = None,
        sample_order: list = None,
        per_sample_overrides: dict = None
    ) -> go.Figure:
        """Create individual graph for each gene with proper data handling"""
        
        # Guard against empty data
        if data is None or data.empty:
            fig = go.Figure()
            fig.add_annotation(text="No data available", showarrow=False)
            return fig
        
        # Work with the gene data directly (it's already filtered by gene)
        gene_data = data.copy()
        
        # Ensure we have required columns
        required_cols = ['Condition', 'Relative_Expression', 'SEM']
        missing = [col for col in required_cols if col not in gene_data.columns]
        
        # Try alternative column names
        if 'Relative_Expression' not in gene_data.columns:
            if 'Fold_Change' in gene_data.columns:
                gene_data['Relative_Expression'] = gene_data['Fold_Change']
            else:
                st.error(f"Missing Relative_Expression or Fold_Change column for {gene}")
                fig = go.Figure()
                fig.add_annotation(text=f"Missing data columns for {gene}", showarrow=False)
                return fig
        
        if 'SEM' not in gene_data.columns:
            gene_data['SEM'] = 0
        
        # ALWAYS use sample_order from mapping (converted to condition names)
        if sample_order:
            # Convert sample names to condition names using mapping
            mapping = st.session_state.get('sample_mapping', {})
            condition_order = []
            for sample in sample_order:
                cond = mapping.get(sample, {}).get('condition', sample)
                if cond in gene_data['Condition'].unique():
                    condition_order.append(cond)
            
            # Add any conditions not in order (shouldn't happen, but safety)
            for cond in gene_data['Condition'].unique():
                if cond not in condition_order:
                    condition_order.append(cond)
            
            # Apply categorical ordering
            gene_data['Condition'] = pd.Categorical(
                gene_data['Condition'], 
                categories=condition_order, 
                ordered=True
            )
            gene_data = gene_data.sort_values('Condition')
        else:
            # Fallback: sort by appearance order
            gene_data = gene_data.sort_values('Condition')
            
        # Get colors - Grey tones for controls, custom for treatments
        # Check for per-sample custom colors first
        bar_colors = []
        
        # Standard grey-tone palette for controls
        control_colors = {
            'Baseline': '#E8E8E8',           # Very light grey (non-treated)
            'Non-treated': '#E8E8E8',        # Very light grey
            'Control': '#E8E8E8',            # Very light grey
            'Negative Control': '#BDBDBD',   # Medium grey
            'Inducer': '#BDBDBD',            # Medium grey (for α-MSH, IGF, etc)
            'Positive Control': '#757575',   # Dark grey
            'Treatment': None                # Will use custom color
        }
        
        for idx, row in gene_data.iterrows():
            condition = row['Condition']
            group = row.get('Group', 'Treatment')
            
            # Check if user has set a custom color for this specific bar
            custom_key = f"{gene}_{condition}"
            if custom_key in settings.get('bar_colors_per_sample', {}):
                bar_colors.append(settings['bar_colors_per_sample'][custom_key])
            # Use control grey tones
            elif group in control_colors and control_colors[group] is not None:
                bar_colors.append(control_colors[group])
            # Use gene default color for treatments
            else:
                default_color = settings.get('bar_colors', {}).get(gene, '#4ECDC4')
                bar_colors.append(default_color)
                
        # Get significance text
        sig_text = gene_data.get('significance', [''] * len(gene_data))
        
        # Create figure
        fig = go.Figure()
        
        # Error bars
        error_array = gene_data['SEM'] * settings.get('error_multiplier', 1.96)
        
        # Add bar trace
        fig.add_trace(go.Bar(
            x=gene_data['Condition'],
            y=gene_data['Relative_Expression'],
            error_y=dict(
                type='data',
                array=error_array,
                visible=settings.get('show_error', True),
                thickness=2,
                width=4,
                color='rgba(0,0,0,0.5)'
            ),
            text=sig_text if settings.get('show_significance', True) else None,
            textposition='outside',
            textfont=dict(size=settings.get('sig_font_size', 16), color='black'),
            marker=dict(
                color=bar_colors,
                line=dict(
                    width=settings.get('marker_line_width', 1),
                    color='black'
                ),
                opacity=settings.get('bar_opacity', 0.95)
            ),
            showlegend=False
        ))
        
        # Add expected direction annotation
        if efficacy_config and 'expected_direction' in efficacy_config:
            direction = efficacy_config.get('expected_direction', {}).get(gene, '')
            if direction == 'up':
                fig.add_annotation(
                    text='↑ Expected Increase',
                    xref='paper', yref='paper',
                    x=0.02, y=0.98,
                    showarrow=False,
                    font=dict(size=12, color='red'),
                    align='left'
                )
            elif direction == 'down':
                fig.add_annotation(
                    text='↓ Expected Decrease',
                    xref='paper', yref='paper',
                    x=0.02, y=0.98,
                    showarrow=False,
                    font=dict(size=12, color='blue'),
                    align='left'
                )
        
        # Update layout
        y_axis_config = dict(
            title=settings.get('ylabel', 'Fold Change'),
            showgrid=settings.get('show_grid', True),
            gridcolor='lightgray'
        )
        
        if settings.get('y_log_scale'):
            y_axis_config['type'] = 'log'
        
        if settings.get('y_min') is not None or settings.get('y_max') is not None:
            y_range = []
            if settings.get('y_min') is not None:
                y_range.append(settings['y_min'])
            if settings.get('y_max') is not None:
                y_range.append(settings['y_max'])
            if len(y_range) == 2:
                y_axis_config['range'] = y_range
        
        # Get gene-specific settings (fallback to global)
        gene_bar_gap = settings.get(f"{gene}_bar_gap", settings.get('bar_gap', 0.15))
        gene_margins = settings.get(f"{gene}_margins", {'l': 80, 'r': 80, 't': 100, 'b': 100})
        gene_bg_color = settings.get(f"{gene}_bg_color", settings.get('plot_bgcolor', '#FFFFFF'))
        gene_grid_color = settings.get(f"{gene}_grid_color", '#E5E5E5')
        gene_grid_width = settings.get(f"{gene}_grid_width", 1)
        gene_xlabel_size = settings.get(f"{gene}_xlabel_size", 14)
        gene_ylabel_size = settings.get(f"{gene}_ylabel_size", 14)
        gene_xlabel_color = settings.get(f"{gene}_xlabel_color", '#000000')
        gene_ylabel_color = settings.get(f"{gene}_ylabel_color", '#000000')
        gene_tick_size = settings.get(f"{gene}_tick_size", 12)
        gene_tick_angle = settings.get(f"{gene}_tick_angle", -45)
        
        fig.update_layout(
            title=dict(
                text=f"{gene} Expression",
                font=dict(size=settings.get('title_size', 20))
            ),
            xaxis=dict(
                title=dict(
                    text=settings.get('xlabel', 'Condition'),
                    font=dict(size=gene_xlabel_size, color=gene_xlabel_color)
                ),
                showgrid=settings.get('show_grid', True),
                gridcolor=gene_grid_color,
                gridwidth=gene_grid_width,
                tickangle=gene_tick_angle,
                tickfont=dict(size=gene_tick_size)
            ),
            yaxis=y_axis_config,
            template=settings.get('color_scheme', 'plotly_white'),
            font=dict(size=settings.get('font_size', 14)),
            height=settings.get('figure_height', 600),
            width=settings.get('figure_width', 1000),
            bargap=gene_bar_gap,
            showlegend=settings.get('show_legend', False),
            plot_bgcolor=gene_bg_color,
            margin=dict(
                l=gene_margins['l'],
                r=gene_margins['r'],
                t=gene_margins['t'],
                b=gene_margins['b']
            )
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
        
        # Raw data (include mapped Condition column reflecting sample mapping)
        raw_export = raw_data.copy()
        # mapping provided as argument 'mapping'
        if mapping:
            raw_export['Condition'] = raw_export['Sample'].map(lambda x: mapping.get(x, {}).get('condition', x))
        else:
            raw_export['Condition'] = raw_export['Sample']
        # Keep original Sample name but add mapped Condition
        raw_export = raw_export[['Well','Sample','Condition','Target','CT','Source_File']] if 'Source_File' in raw_export.columns else raw_export
        raw_export.to_excel(writer, sheet_name='Raw_Data', index=False)

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
        
        # Sample mapping interface with professional layout
        st.markdown("---")
        st.markdown("### 🗺️ Sample Condition Mapping")
        
        samples = [s for s in sorted(st.session_state.data['Sample'].unique()) 
                  if s not in st.session_state.excluded_samples]
        
        # Group type options
        group_types = ['Negative Control', 'Positive Control', 'Treatment']
        if 'baseline' in config['controls']:
            group_types.insert(0, 'Baseline')
        
        # Initialize sample_order if not exists
        if 'sample_order' not in st.session_state:
            st.session_state.sample_order = samples.copy()
        
        # Ensure all samples have mapping
        for sample in samples:
            if sample not in st.session_state.sample_mapping:
                st.session_state.sample_mapping[sample] = {
                    'condition': sample, 
                    'group': 'Treatment', 
                    'concentration': '', 
                    'include': True
                }
            if 'include' not in st.session_state.sample_mapping[sample]:
                st.session_state.sample_mapping[sample]['include'] = True
        
        # Master controls in styled container
        st.markdown("""
        <style>
        .control-panel {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            padding: 20px;
            border-radius: 10px;
            margin-bottom: 20px;
        }
        .control-button {
            background-color: white;
            color: #667eea;
            border: none;
            padding: 10px 20px;
            border-radius: 5px;
            cursor: pointer;
            font-weight: bold;
        }
        </style>
        """, unsafe_allow_html=True)
        
        st.markdown('<div class="control-panel">', unsafe_allow_html=True)
        col_master1, col_master2, col_master3 = st.columns(3)
        
        with col_master1:
            if st.button("✅ Include ALL Samples", use_container_width=True):
                for s in st.session_state.sample_order:
                    st.session_state.sample_mapping[s]['include'] = True
                st.session_state.excluded_samples = set()
                st.rerun()
        
        with col_master2:
            if st.button("🚫 Exclude ALL Samples", use_container_width=True):
                for s in st.session_state.sample_order:
                    st.session_state.sample_mapping[s]['include'] = False
                st.session_state.excluded_samples = set(st.session_state.sample_order)
                st.rerun()
        
        with col_master3:
            if st.button("🔄 Reset Mapping", use_container_width=True):
                st.session_state.sample_mapping = {}
                st.session_state.sample_order = samples.copy()
                st.rerun()
        
        st.markdown('</div>', unsafe_allow_html=True)
        
        # Header row with styled background
        st.markdown("""
        <div style='background-color: #f8f9fa; padding: 10px; border-radius: 5px; margin-bottom: 10px;'>
            <table style='width: 100%;'>
                <tr>
                    <th style='width: 5%; text-align: center;'>✓</th>
                    <th style='width: 10%;'>Order</th>
                    <th style='width: 15%;'>Original</th>
                    <th style='width: 25%;'>Condition Name</th>
                    <th style='width: 20%;'>Group</th>
                    <th style='width: 15%;'>Concentration</th>
                    <th style='width: 10%;'>Move</th>
                </tr>
            </table>
        </div>
        """, unsafe_allow_html=True)
        
        # Sample rows with improved spacing
        for i, sample in enumerate(st.session_state.sample_order):
            # Container for each row
            with st.container():
                col0, col_order, col1, col2, col3, col4, col_move = st.columns([0.5, 0.8, 1.5, 2.5, 2, 1.5, 1])
                
                # Include checkbox
                with col0:
                    include = st.checkbox(
                        "", 
                        value=st.session_state.sample_mapping[sample].get('include', True),
                        key=f"include_{sample}",
                        label_visibility="collapsed"
                    )
                    st.session_state.sample_mapping[sample]['include'] = include
                
                # Order number
                with col_order:
                    st.markdown(f"<div style='text-align: center; padding-top: 10px;'><b>{i+1}</b></div>", unsafe_allow_html=True)
                
                # Original sample name (non-editable)
                with col1:
                    st.text_input("Original", sample, key=f"orig_{sample}", disabled=True, label_visibility="collapsed")
                
                # Condition name (editable)
                with col2:
                    cond = st.text_input(
                        "Condition",
                        st.session_state.sample_mapping[sample]['condition'],
                        key=f"cond_{sample}",
                        label_visibility="collapsed",
                        placeholder="Enter condition name..."
                    )
                    st.session_state.sample_mapping[sample]['condition'] = cond
                
                # Group selector
                with col3:
                    grp_idx = 0
                    try:
                        grp_idx = group_types.index(st.session_state.sample_mapping[sample]['group'])
                    except:
                        pass
                    
                    grp = st.selectbox(
                        "Group",
                        group_types,
                        index=grp_idx,
                        key=f"grp_{sample}",
                        label_visibility="collapsed"
                    )
                    st.session_state.sample_mapping[sample]['group'] = grp
                
                # Concentration
                with col4:
                    conc = st.text_input(
                        "Conc",
                        st.session_state.sample_mapping[sample]['concentration'],
                        key=f"conc_{sample}",
                        label_visibility="collapsed",
                        placeholder="e.g., 10µM"
                    )
                    st.session_state.sample_mapping[sample]['concentration'] = conc
                
                # Move buttons
                with col_move:
                    move_col1, move_col2 = st.columns(2)
                    with move_col1:
                        if i > 0:
                            if st.button("⬆", key=f"up_{sample}", help="Move up"):
                                order = st.session_state.sample_order
                                order[i-1], order[i] = order[i], order[i-1]
                                st.rerun()
                    with move_col2:
                        if i < len(st.session_state.sample_order) - 1:
                            if st.button("⬇", key=f"down_{sample}", help="Move down"):
                                order = st.session_state.sample_order
                                order[i+1], order[i] = order[i], order[i+1]
                                st.rerun()
                
                # Divider line
                st.markdown("<hr style='margin: 5px 0; opacity: 0.3;'>", unsafe_allow_html=True)
        
        # Update excluded_samples from include flags
        st.session_state.excluded_samples = set([
            s for s, v in st.session_state.sample_mapping.items() 
            if not v.get('include', True)
        ])
        
        # Summary with styled cards
        st.markdown("---")
        st.subheader("📊 Mapping Summary")
        
        col_card1, col_card2, col_card3, col_card4 = st.columns(4)
        
        total_samples = len(st.session_state.sample_order)
        included = sum(1 for v in st.session_state.sample_mapping.values() if v.get('include', True))
        excluded = total_samples - included
        
        with col_card1:
            st.metric("Total Samples", total_samples)
        with col_card2:
            st.metric("Included", included, delta=None if included == total_samples else f"-{excluded}")
        with col_card3:
            st.metric("Excluded", excluded, delta=None if excluded == 0 else f"+{excluded}")
        with col_card4:
            groups = set(v['group'] for v in st.session_state.sample_mapping.values() if v.get('include', True))
            st.metric("Groups", len(groups))
        
        # Detailed table view
        with st.expander("📋 View Detailed Mapping Table"):
            mapping_df = pd.DataFrame([
                {
                    'Order': idx+1,
                    'Include': '✅' if v.get('include', True) else '❌',
                    'Original': s,
                    'Condition': v['condition'],
                    'Group': v['group'],
                    'Concentration': v['concentration']
                }
                for idx, s in enumerate(st.session_state.sample_order)
                for v in [st.session_state.sample_mapping[s]]
            ])
            st.dataframe(mapping_df, use_container_width=True, hide_index=True)
          
        # Run analysis   
        st.markdown("---")
        st.subheader("🔬 Run Full Analysis (ΔΔCt + Statistics)")
        
        sample_keys = st.session_state.get('sample_order') or list(st.session_state.sample_mapping.keys())
        
        if sample_keys:
            # Enhanced layout with clear separation
            st.markdown("#### 📊 Analysis Configuration")
            
            col_info1, col_info2 = st.columns(2)
            with col_info1:
                st.info("**ΔΔCt Reference:** Used to calculate fold changes. All samples will be relative to this (Fold Change = 1.0)")
            with col_info2:
                st.info("**P-value Reference:** Used for statistical comparison (t-test). All samples compared against this.")
            
            col_r1, col_r2 = st.columns(2)
            with col_r1:
                ref_choice = st.selectbox(
                    "🎯 ΔΔCt Reference Sample",
                    sample_keys,
                    index=0,
                    key="ref_choice_ddct",
                    help="Baseline for relative expression calculation"
                )
                ref_condition = st.session_state.sample_mapping.get(ref_choice, {}).get('condition', ref_choice)
                st.caption(f"→ Condition: **{ref_condition}**")
            
            with col_r2:
                cmp_choice = st.selectbox(
                    "📈 P-value Reference Sample",
                    sample_keys,
                    index=0,
                    key="cmp_choice_pval",
                    help="Control group for statistical testing"
                )
                cmp_condition = st.session_state.sample_mapping.get(cmp_choice, {}).get('condition', cmp_choice)
                st.caption(f"→ Condition: **{cmp_condition}**")
            
            # Visual summary
            st.markdown("---")
            col_sum1, col_sum2, col_sum3 = st.columns([1, 2, 1])
            with col_sum2:
                st.markdown(f"""
                <div style='background-color: #f0f2f6; padding: 15px; border-radius: 10px; text-align: center;'>
                    <h4>Analysis Summary</h4>
                    <p><b>Fold Changes:</b> Relative to <code>{ref_condition}</code></p>
                    <p><b>Significance:</b> Compared to <code>{cmp_condition}</code></p>
                </div>
                """, unsafe_allow_html=True)
            
            st.markdown("---")
            
            # Run button
            if st.button("▶️ Run Full Analysis Now", type="primary", use_container_width=True):
                ok = AnalysisEngine.run_full_analysis(ref_choice, cmp_choice)
                if ok:
                    st.success(f"✅ Analysis complete!\n\n- Fold changes relative to: **{ref_condition}**\n- P-values vs: **{cmp_condition}**")
                    st.balloons()
                    st.rerun()
                else:
                    st.error("❌ Analysis failed. Check messages above.")
        else:
            st.warning("⚠️ No samples available for analysis.")
            
# ==================== TAB 3: ANALYSIS ====================
with tab3:
    st.header("Step 3: Analysis Results")
    
    if st.session_state.processed_data:
        st.subheader("📊 Analysis Summary")
        
        # Summary metrics
        all_results = pd.concat(st.session_state.processed_data.values(), ignore_index=True)
        
        col1, col2, col3 = st.columns(3)
        col1.metric("Genes Analyzed", len(st.session_state.processed_data))
        col2.metric("Conditions", all_results['Condition'].nunique())
        sig_count = (all_results['p_value'] < 0.05).sum()
        col3.metric("Significant (p<0.05)", f"{sig_count}/{len(all_results)}")
        
        # Show results per gene
        st.subheader("🧬 Gene-by-Gene Results")
        
        for gene, gene_df in st.session_state.processed_data.items():
            with st.expander(f"📍 {gene}", expanded=False):
                # Show expected direction if available
                efficacy_config = EFFICACY_CONFIG.get(st.session_state.selected_efficacy, {})
                if 'expected_direction' in efficacy_config:
                    direction = efficacy_config['expected_direction'].get(gene)
                    if direction:
                        st.caption(f"Expected: {'↑ Increase' if direction == 'up' else '↓ Decrease'}")
                
                # Display columns
                display_cols = ['Condition', 'Group', 'Fold_Change', 'p_value', 'significance', 
                              'n_replicates', 'Target_Ct_Mean', 'HK_Ct_Mean', 'Delta_Ct', 'SEM']
                
                # Filter to existing columns
                display_df = gene_df[[c for c in display_cols if c in gene_df.columns]]
                
                # Style the dataframe
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
        
        st.success("✅ Results ready! Go to Graphs tab to visualize.")
    
    else:
        st.info("⏳ No analysis results yet. Go to 'Sample Mapping' tab and click 'Run Full Analysis Now'")
        
# ==================== TAB 4: GRAPHS (IMPROVED UI/UX - ERROR-SAFE) ====================
with tab4:
    st.header("Step 4: Advanced Gene Visualization 📈")
    st.markdown("---")
    
    if st.session_state.processed_data:
        
        # --- Persistent Settings Initialization ---
        # 1. Global Settings
        if 'graph_settings' not in st.session_state:
            st.session_state.graph_settings = { 
                'title_size': 24, 'font_size': 12, 'sig_font_size': 18, 
                'figure_width': 700, 'figure_height': 500, 
                'color_scheme': 'plotly_white', 'show_error': True, 
                'show_significance': True, 'show_grid': True, 
                'xlabel': 'Condition', 'ylabel': 'Fold Change (Relative Expression)', 
                'bar_colors': {}, 'error_multiplier': 1.96, 'bar_opacity': 0.95, 
                'bar_gap': 0.1, 'marker_line_width': 1, 'show_legend': False, 
                'y_log_scale': False, 'y_min': 0, 'y_max': None,
                'bar_colors_per_sample': {}
            }
            
        # 2. Per-Gene Customizations (Dedicated Dictionary for safety and clarity)
        if 'gene_customizations' not in st.session_state:
            st.session_state.gene_customizations = {}

        # --- PRE-INITIALIZE ALL GENE SETTINGS (CRITICAL ERROR FIX) ---
        # This loop runs once or on rerun, ensuring all necessary keys exist 
        # BEFORE the main plotting loop (the one that causes the error) begins.
        for gene in st.session_state.processed_data.keys():
            if gene not in st.session_state.gene_customizations:
                st.session_state.gene_customizations[gene] = {
                    'title_size': st.session_state.graph_settings['title_size'],
                    'y_log_scale': st.session_state.graph_settings['y_log_scale'],
                    'y_min': st.session_state.graph_settings['y_min'],
                    'y_max': st.session_state.graph_settings['y_max'],
                    'bar_gap': st.session_state.graph_settings['bar_gap'],
                    'tick_angle': -45, # Default value
                    'bg_color': '#FFFFFF', # Default value
                    'grid_color': '#E5E5E5', # Default value
                }
        
        # --- Global Defaults/Reset Panel (Compact Expander) ---
        # This replaces the old, large Global Graph Settings panel
        with st.expander("🛠️ Global Visualization Defaults & Reset", expanded=False):
            col_g1, col_g2, col_g3 = st.columns(3)
            
            with col_g1:
                st.session_state.graph_settings['figure_width'] = st.number_input(
                    "Default Width (px)", 300, 2000, st.session_state.graph_settings['figure_width'], key='g_width', step=100
                )
                st.session_state.graph_settings['figure_height'] = st.number_input(
                    "Default Height (px)", 300, 1400, st.session_state.graph_settings['figure_height'], key='g_height', step=100
                )
            
            with col_g2:
                themes = ['plotly_white', 'plotly', 'simple_white', 'plotly_dark', 'seaborn', 'presentation']
                curr_theme = st.session_state.graph_settings.get('color_scheme', 'plotly_white')
                st.session_state.graph_settings['color_scheme'] = st.selectbox(
                    "Default Theme", themes, index=themes.index(curr_theme) if curr_theme in themes else 0, key='g_theme'
                )
                st.session_state.graph_settings['font_size'] = st.slider( "Default Font Size", 8, 20, st.session_state.graph_settings['font_size'], key='g_font' )
            
            with col_g3:
                st.session_state.graph_settings['show_error'] = st.checkbox("Show Error Bars", st.session_state.graph_settings['show_error'], key='g_error')
                st.session_state.graph_settings['show_significance'] = st.checkbox("Show Significance Stars", st.session_state.graph_settings['show_significance'], key='g_sig')
                if st.button("🔄 Reset ALL Settings", use_container_width=True):
                    # Safely delete keys and force a rerun
                    for key in ['graph_settings', 'gene_customizations', 'graphs']:
                         if key in st.session_state:
                             del st.session_state[key]
                    st.rerun()

        st.markdown("---")
        
        # --- Gene-by-Gene Iteration (New UI/UX) ---
        efficacy_config = EFFICACY_CONFIG.get(st.session_state.selected_efficacy, {})
        
        # Use columns for a clean 2-column layout of genes
        gene_columns = st.columns(2)
        
        # Iterate over the stable dictionary keys
        for i, gene in enumerate(st.session_state.processed_data.keys()):
            # Get the dedicated settings dictionary for the current gene
            gene_settings = st.session_state.gene_customizations[gene] 
            
            with gene_columns[i % 2]:
                
                # --- High-Tech Panel Container (CSS for WoW factor) ---
                st.markdown(f"""
                <div style='
                    background-color: #f8f9fa; 
                    padding: 15px; 
                    border-radius: 10px; 
                    margin-bottom: 20px; 
                    box-shadow: 0 4px 8px rgba(0,0,0,0.1);
                    border-left: 5px solid #007bff;
                '>
                    <h3 style='margin-top: 0; color: #007bff;'>🧬 {gene} Result</h3>
                """, unsafe_allow_html=True)
                
                # --- 1. Compact Customization Bar ---
                col_c1, col_c2, col_c3, col_c4, col_c5 = st.columns([1, 1.5, 1, 1, 1])
                
                # Toggle Log Scale (Uses gene_settings directly)
                with col_c1:
                    gene_settings['y_log_scale'] = st.checkbox(
                        "Log Y", gene_settings['y_log_scale'], key=f"log_y_{gene}", help="Toggle Y-axis logarithmic scale", label_visibility="collapsed"
                    )
                    st.caption("Log Y")
                    
                # Title Size Override (Uses gene_settings directly)
                with col_c2:
                    gene_settings['title_size'] = st.slider(
                        "Title Size", 10, 40, gene_settings['title_size'], key=f"ts_{gene}", label_visibility="collapsed"
                    )
                    st.caption("Title Size")
                
                # Y-Range Popover
                with col_c3:
                    with st.popover("📐 Range"):
                        st.markdown("**Y-Axis Range**")
                        
                        gene_settings['y_min'] = st.number_input(
                            "Min Y", value=gene_settings['y_min'], key=f"ymin_{gene}"
                        )
                        gene_settings['y_max'] = st.number_input(
                            "Max Y", value=gene_settings['y_max'], key=f"ymax_{gene}"
                        )
                    st.caption("Y-Axis")

                # Advanced Layout Popover
                with col_c4:
                    with st.popover("⚙️ Layout"):
                        st.markdown("**Advanced Layout**")
                        col_l1, col_l2 = st.columns(2)
                        
                        with col_l1:
                            # Bar Gap
                            gene_settings['bar_gap'] = st.slider(
                                "Bar Gap", 0.0, 0.5, gene_settings['bar_gap'], 0.01, key=f"gap_{gene}", label_visibility="collapsed"
                            )
                            # X-Tick Angle
                            gene_settings['tick_angle'] = st.slider( 
                                "X-Tick Angle", -90, 90, gene_settings['tick_angle'], key=f"tickang_{gene}", label_visibility="collapsed"
                            )
                            
                        with col_l2:
                            # Plot BG Color
                            gene_settings['bg_color'] = st.color_picker(
                                "BG Color", gene_settings['bg_color'], key=f"bgc_{gene}", label_visibility="collapsed"
                            )
                            # Grid Color
                            gene_settings['grid_color'] = st.color_picker(
                                "Grid Color", gene_settings['grid_color'], key=f"gc_{gene}", label_visibility="collapsed"
                            )
                    st.caption("Layout")
                
                # Color Picker Popover
                with col_c5:
                    gene_data_for_color = st.session_state.processed_data[gene]
                    conditions = gene_data_for_color['Condition'].unique()
                    
                    with st.popover("🎨 Colors"):
                        st.markdown(f"**Gene Default Color:**")
                        default_color = st.session_state.graph_settings['bar_colors'].get(gene, '#007bff')
                        st.session_state.graph_settings['bar_colors'][gene] = st.color_picker(
                            f"Default Color for {gene}", default_color, key=f"gene_color_{gene}"
                        )
                        st.markdown("---")
                        st.markdown(f"**Individual Bar Overrides:**")
                        
                        # Loop through conditions to offer individual overrides
                        for condition in conditions:
                            custom_key = f"{gene}_{condition}"
                            
                            # Determine current color (default to gene default if no override)
                            current_bar_color = st.session_state.graph_settings['bar_colors_per_sample'].get(
                                custom_key, 
                                st.session_state.graph_settings['bar_colors'][gene]
                            )
                            
                            col_bar_c, col_bar_r = st.columns([3, 1])
                            
                            with col_bar_c:
                                new_color = st.color_picker( 
                                    f"Bar: {condition}", current_bar_color, key=f"bar_color_{gene}_{condition}" 
                                )
                                # Save override color
                                st.session_state.graph_settings['bar_colors_per_sample'][custom_key] = new_color
                            
                            with col_bar_r:
                                st.markdown("<div style='padding-top: 25px;'>", unsafe_allow_html=True)
                                # The reset button deletes from a nested dictionary, which is safe.
                                if st.button("X", key=f"reset_{gene}_{condition}", help="Reset to gene default", use_container_width=True):
                                    if custom_key in st.session_state.graph_settings['bar_colors_per_sample']:
                                        del st.session_state.graph_settings['bar_colors_per_sample'][custom_key]
                                    st.rerun() 

                    st.caption("Bar Colors")


                st.markdown("<hr style='border: 1px dashed #ced4da;'>", unsafe_allow_html=True)

                # --- 2. Generate and Display Graph ---
                
                # Merge global settings with gene-specific overrides for final plotting
                final_settings = st.session_state.graph_settings.copy()
                final_settings.update(gene_settings) # Override global settings with gene-specific settings
                
                gene_data = st.session_state.processed_data[gene]
                sample_order = st.session_state.get('sample_order')
                
                fig = GraphGenerator.create_gene_graph(
                    gene_data,
                    gene,
                    final_settings, # Pass the consolidated final settings
                    efficacy_config,
                    sample_order
                )
                
                # Final touch: store the generated figure for export
                st.session_state.graphs[gene] = fig
                st.plotly_chart(fig, use_container_width=True)
                
                # --- Panel Footer (Close the container) ---
                st.markdown("</div>", unsafe_allow_html=True)

    else:
        st.warning("⚠️ Run the analysis in the 'Analysis' tab first.")
    
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
            'Reference_Sample': 'ref_condition' if 'ref_condition' in locals() else 'N/A',
            'Compare_To': 'compare_condition' if 'compare_condition' in locals() else 'N/A',
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

