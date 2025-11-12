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
    def calculate_statistics(processed: pd.DataFrame, compare_condition: str,
                             raw_data: pd.DataFrame = None,
                             hk_gene: str = None,
                             sample_mapping: dict = None) -> pd.DataFrame:
        """
        Perform two-tailed Welch's t-test (independent samples) comparing each condition
        to the specified compare_condition. Works safely with triplicates.
        If raw_data, hk_gene, or sample_mapping are not provided, it falls back to
        st.session_state variables. Results are added as p_value and significance columns.
        """
        import numpy as np
        from scipy import stats
        import streamlit as st

        # Use session_state fallbacks if args not provided
        raw_data = raw_data or st.session_state.get("data")
        hk_gene = hk_gene or st.session_state.get("hk_gene")
        sample_mapping = sample_mapping or st.session_state.get("sample_mapping", {})

        # If raw_data missing, return processed unchanged
        if raw_data is None or hk_gene is None:
            return processed

        results = processed.copy()
        results["p_value"] = np.nan
        results["significance"] = ""

        for target in results["Target"].unique():
            if pd.isna(target):
                continue

            # Map samples to user-defined condition names
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

            rel_expr = {}
            for cond, grp in t_rows.groupby("Condition"):
                hk_mean = hk_means.get(cond, np.nan)
                if np.isnan(hk_mean):
                    continue
                rel_expr[cond] = 2 ** (-(grp["CT"].values - hk_mean))

            # Reference condition replicates
            ref_vals = rel_expr.get(compare_condition, np.array([]))
            if ref_vals.size == 0:
                continue

            # Compare each condition to reference
            for cond, vals in rel_expr.items():
                if cond == compare_condition or vals.size == 0:
                    continue

                try:
                    if len(ref_vals) >= 2 and len(vals) >= 2:
                        _, p_val = stats.ttest_ind(ref_vals, vals, equal_var=False)
                    elif len(ref_vals) == 1 and len(vals) > 1:
                        _, p_val = stats.ttest_1samp(vals, ref_vals[0])
                    elif len(vals) == 1 and len(ref_vals) > 1:
                        _, p_val = stats.ttest_1samp(ref_vals, vals[0])
                    else:
                        p_val = np.nan
                except Exception:
                    p_val = np.nan

                # Annotate results
                mask = (results["Target"] == target) & (results["Condition"] == cond)
                results.loc[mask, "p_value"] = float(p_val) if not np.isnan(p_val) else np.nan
                if not np.isnan(p_val):
                    if p_val < 0.001:
                        sig = "***"
                    elif p_val < 0.01:
                        sig = "**"
                    elif p_val < 0.05:
                        sig = "*"
                    else:
                        sig = ""
                    results.loc[mask, "significance"] = sig

        return results

    def run_full_analysis(ref_sample_key: str, compare_sample_key: str):
        """
        Run a full ΔΔCt + statistics pipeline and store results in
        st.session_state.processed_data as a dict: {gene: DataFrame}.
        - ref_sample_key / compare_sample_key are the ORIGINAL sample keys
        (i.e. keys of st.session_state.sample_mapping), not the mapped condition names.
        """
        try:
            # basic safety checks
            if st.session_state.get('data') is None:
                st.error("No raw data loaded.")
                return False
            if not st.session_state.get('sample_mapping'):
                st.error("No sample mapping available.")
                return False
            if not st.session_state.get('hk_gene'):
                st.error("No housekeeping gene selected.")
                return False

            # Resolve mapped condition names (what calculate_ddct expects in your current file)
            mapping = st.session_state.sample_mapping
            ref_condition = mapping.get(ref_sample_key, {}).get('condition', ref_sample_key)
            compare_condition = mapping.get(compare_sample_key, {}).get('condition', compare_sample_key)

            # Call ΔΔCt calculation (signature in your file: data, hk_gene, ref_sample, compare_sample, excluded_wells, excluded_samples, sample_mapping)
            processed_df = AnalysisEngine.calculate_ddct(
                st.session_state.data,
                st.session_state.hk_gene,
                ref_condition,
                compare_condition,
                st.session_state.get('excluded_wells', set()),
                st.session_state.get('excluded_samples', set()),
                mapping
            )

            if processed_df is None or processed_df.empty:
                st.warning("No processed results (ΔΔCt) — please check inputs.")
                return False

            # Calculate statistics. Try to call the robust signature but fall back safely.
            try:
                processed_with_stats = AnalysisEngine.calculate_statistics(
                    processed_df,
                    compare_condition,
                    raw_data=st.session_state.data,
                    hk_gene=st.session_state.hk_gene,
                    sample_mapping=mapping
                )
            except TypeError:
                # older signature: calculate_statistics(processed, compare_condition)
                processed_with_stats = AnalysisEngine.calculate_statistics(processed_df, compare_condition)

            # Convert to dict of dataframes keyed by gene (this is what the graphs tab expects)
            gene_dict = {}
            if 'Target' in processed_with_stats.columns:
                for g in processed_with_stats['Target'].unique():
                    gene_df = processed_with_stats[processed_with_stats['Target'] == g].copy()
                    gene_dict[g] = gene_df.reset_index(drop=True)
            else:
                # If output already a dict-like, try to use it directly (defensive)
                if isinstance(processed_with_stats, dict):
                    gene_dict = processed_with_stats
                else:
                    # fallback: single table under a generic key
                    gene_dict = {'results': processed_with_stats.reset_index(drop=True)}

            st.session_state.processed_data = gene_dict
            return True

        except Exception as e:
            st.error(f"Analysis failed: {e}")
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
        """
        Create an interactive Plotly bar graph for the given gene, using user-defined settings.
        Handles empty data safely and supports sample ordering + per-sample overrides.
        """

        # Guard against empty data
        if data is None or data.empty:
            fig = go.Figure()
            fig.add_annotation(text="No data available for plotting", showarrow=False)
            return fig

        # Filter for this gene
        if "Target" in data.columns:
            gene_data = data[data["Target"] == gene].copy()
        elif "Gene" in data.columns:
            gene_data = data[data["Gene"] == gene].copy()
        else:
            gene_data = data.copy()

        if gene_data.empty:
            fig = go.Figure()
            fig.add_annotation(text=f"No data found for {gene}", showarrow=False)
            return fig

        # Maintain custom sample order
        if sample_order:
            gene_data["Condition"] = pd.Categorical(
                gene_data["Condition"], categories=sample_order, ordered=True
            )
            gene_data = gene_data.sort_values("Condition")

        # Apply per-sample visibility overrides
        if per_sample_overrides:
            keep_rows = []
            for _, row in gene_data.iterrows():
                key = str((gene, row["Condition"]))
                include = per_sample_overrides.get(key, {}).get("include", True)
                keep_rows.append(include)
            gene_data = gene_data[np.array(keep_rows)]

        # Apply bar colors
        default_color = settings.get("bar_colors", {}).get(
            gene, "#4ECDC4"
        )
        bar_colors = []
        for _, row in gene_data.iterrows():
            key = str((gene, row["Condition"]))
            color = (
                per_sample_overrides.get(key, {}).get("color")
                if per_sample_overrides
                else None
            )
            bar_colors.append(color or default_color)

        # Ensure numeric values
        y_vals = pd.to_numeric(gene_data["Relative_Expression"], errors="coerce")
        err_vals = pd.to_numeric(gene_data.get("SEM", 0), errors="coerce")

        # Build figure
        fig = go.Figure(
            data=[
                go.Bar(
                    x=gene_data["Condition"],
                    y=y_vals,
                    error_y=dict(
                        type="data",
                        array=err_vals,
                        visible=settings.get("show_error", True),
                    ),
                    text=gene_data.get("significance", ""),
                    textposition="outside",
                    marker=dict(
                        color=bar_colors,
                        line=dict(width=1, color="black"),
                    ),
                )
            ]
        )

        # Annotate expected direction if configured
        if efficacy_config and "expected_direction" in efficacy_config:
            direction = efficacy_config["expected_direction"].get(gene)
            if direction == "up":
                fig.add_annotation(
                    text="↑ Expected Increase",
                    xref="paper", yref="paper",
                    x=0.02, y=0.95,
                    showarrow=False, font=dict(color="red", size=12)
                )
            elif direction == "down":
                fig.add_annotation(
                    text="↓ Expected Decrease",
                    xref="paper", yref="paper",
                    x=0.02, y=0.95,
                    showarrow=False, font=dict(color="blue", size=12)
                )

        # Layout settings
        fig.update_layout(
            title=f"{gene} Expression",
            xaxis_title=settings.get("xlabel", "Condition"),
            yaxis_title=settings.get("ylabel", "Fold Change"),
            template=settings.get("color_scheme", "plotly_white"),
            font=dict(size=settings.get("font_size", 14)),
            height=settings.get("figure_height", 600),
            width=settings.get("figure_width", 1000),
            yaxis=dict(
                type="log" if settings.get("y_log_scale", False) else "linear",
                range=[
                    settings.get("y_min", None),
                    settings.get("y_max", None),
                ],
            ),
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
        
        # Sample mapping interface
        st.subheader("🗺️ Sample Mapping")
        
        samples = [s for s in sorted(st.session_state.data['Sample'].unique()) 
                  if s not in st.session_state.excluded_samples]
        
        # Group type options based on efficacy
        group_types = ['Negative Control', 'Positive Control', 'Treatment']
        if 'baseline' in config['controls']:
            group_types.insert(0, 'Baseline')
        
        # ---- Improved Sample Mapping UI: include/exclude + ordering ----
        # Initialize persistent list for ordering
        if 'sample_order' not in st.session_state:
            st.session_state.sample_order = [s for s in sorted(st.session_state.data['Sample'].unique()) 
                                            if s not in st.session_state.excluded_samples]

        # Master include/exclude buttons
        col_a, col_b, col_c = st.columns([1,1,2])
        with col_a:
            if st.button("✅ Include ALL"):
                for s in st.session_state.sample_order:
                    st.session_state.sample_mapping.setdefault(s, {})
                    st.session_state.sample_mapping[s]['include'] = True
                st.session_state.excluded_samples = set()
                st.rerun()
        with col_b:
            if st.button("🚫 Exclude ALL"):
                for s in st.session_state.sample_order:
                    st.session_state.sample_mapping.setdefault(s, {})
                    st.session_state.sample_mapping[s]['include'] = False
                st.session_state.excluded_samples = set(st.session_state.sample_order)
                st.rerun()
        with col_c:
            st.markdown("Use the move buttons to reorder samples; that order is used for plotting.")

        # Ensure every sample has mapping keys and include flag
        for sample in list(st.session_state.sample_order):
            if sample not in st.session_state.sample_mapping:
                st.session_state.sample_mapping[sample] = {'condition': sample, 'group': 'Treatment', 'concentration': '', 'include': True}
            if 'include' not in st.session_state.sample_mapping[sample]:
                st.session_state.sample_mapping[sample]['include'] = True

        # Render each sample with include checkbox, mapping fields and order controls
        for i, sample in enumerate(st.session_state.sample_order):
            col0, col1, col2, col3, col4 = st.columns([0.6, 1.5, 3, 2, 1.2])
            # include toggle
            include_flag = col0.checkbox("", value=st.session_state.sample_mapping[sample].get('include', True), key=f"include_{sample}")
            st.session_state.sample_mapping[sample]['include'] = include_flag
            # show sample name
            col1.text(sample)
            # condition name
            cond = col2.text_input("Condition name", value=st.session_state.sample_mapping[sample].get('condition', sample), key=f"cond_{sample}", label_visibility="collapsed")
            st.session_state.sample_mapping[sample]['condition'] = cond
            # group selector
            grp_idx = 0
            try:
                grp_idx = group_types.index(st.session_state.sample_mapping[sample].get('group', 'Treatment'))
            except Exception:
                grp_idx = 0
            grp = col3.selectbox("Group", group_types, index=grp_idx, key=f"grp_{sample}", label_visibility="collapsed")
            st.session_state.sample_mapping[sample]['group'] = grp
            # conc
            conc = col4.text_input("Conc", value=st.session_state.sample_mapping[sample].get('concentration',''), key=f"conc_{sample}", label_visibility="collapsed")
            st.session_state.sample_mapping[sample]['concentration'] = conc

            # Order controls (Move Up / Move Down) placed after the row for readability
            col_left, col_right = st.columns([1,1])
            with col_left:
                if st.button("⬆ Move Up", key=f"up_{sample}") and i>0:
                    order = st.session_state.sample_order
                    order[i-1], order[i] = order[i], order[i-1]
                    st.session_state.sample_order = order
                    st.rerun()
            with col_right:
                if st.button("⬇ Move Down", key=f"down_{sample}") and i < len(st.session_state.sample_order)-1:
                    order = st.session_state.sample_order
                    order[i+1], order[i] = order[i], order[i+1]
                    st.session_state.sample_order = order
                    st.rerun()

        # Update excluded_samples set from per-sample include flags
        st.session_state.excluded_samples = set([s for s,v in st.session_state.sample_mapping.items() if not v.get('include', True)])

        # Summary (now includes order)
        st.subheader("📊 Mapping Summary (ordered)")
        mapping_df = pd.DataFrame([{'Order': idx+1, 'Original': k, **v} 
                                for idx, (k,v) in enumerate([(s, st.session_state.sample_mapping[s]) for s in st.session_state.sample_order])])
        st.dataframe(mapping_df, use_container_width=True)
        
        # --- Analysis quick-run controls (placed after mapping summary) ---
        st.markdown("---")
        st.subheader("▶ Run Full Analysis (ΔΔCt + statistics)")

        # Build choices from sample_mapping keys (preserve user ordering if available)
        sample_keys = st.session_state.get('sample_order') or list(st.session_state.sample_mapping.keys())
        if not sample_keys:
            sample_keys = list(st.session_state.sample_mapping.keys())

        col_r1, col_r2 = st.columns(2)
        with col_r1:
            ref_choice = st.selectbox(
                "Reference sample (normalization ΔΔCt)",
                sample_keys,
                index=0 if sample_keys else 0,
                help="This mapped sample will be used as the ΔΔCt baseline (fold change = 1.0)."
            )
        with col_r2:
            cmp_choice = st.selectbox(
                "Comparison sample (for p-values)",
                sample_keys,
                index=0 if sample_keys else 0,
                help="This sample's replicates are used as the comparison group when calculating p-values."
            )

        run_col1, run_col2 = st.columns([1, 3])
        with run_col1:
            if st.button("🔬 Run Full Analysis Now", type="primary"):
                ok = 'run_full_analysis'(ref_choice, cmp_choice)
                if ok:
                    st.success("Analysis finished — go to the Graphs tab to visualize results.")
                    st.rerun()
                else:
                    st.error("Analysis did not complete successfully. Check logs above.")
        with run_col2:
            st.markdown(
                "Tip: After mapping is done, use 'Run Full Analysis Now' to compute ΔΔCt & p-values. "
                "Results will appear in the Graphs and Export tabs."
            )


# ==================== TAB 3: ANALYSIS ====================
with tab3:
    st.header("📊 Visualization")

    # Guard: Require processed data
    if "processed_data" not in st.session_state or st.session_state.processed_data is None:
        st.info("Run analysis first to visualize results.")
    else:
        st.divider()
        st.markdown("### ⚙️ Global Display Settings")

        # Initialize graph settings if missing
        if "graph_settings" not in st.session_state:
            st.session_state.graph_settings = {
                "bar_colors": {},
                "color_scheme": "plotly_white",
                "show_error": True,
                "show_significance": True,
                "y_log_scale": False,
                "show_grid": True,
                "bar_gap": 0.15,
                "font_size": 14,
                "figure_height": 600,
                "figure_width": 900
            }

        # --- Visualization Controls ---
        with st.expander("🎨 Customize Global Graph Appearance", expanded=False):
            c1, c2, c3 = st.columns(3)
            with c1:
                st.session_state.graph_settings["color_scheme"] = st.selectbox(
                    "Color Theme",
                    ["plotly_white", "plotly_dark", "ggplot2", "seaborn", "simple_white"],
                    index=0,
                )
                st.session_state.graph_settings["show_error"] = st.checkbox(
                    "Show Error Bars", value=True
                )
                st.session_state.graph_settings["show_significance"] = st.checkbox(
                    "Show Significance", value=True
                )
            with c2:
                st.session_state.graph_settings["y_log_scale"] = st.checkbox(
                    "Use Log Scale (Y)", value=False
                )
                st.session_state.graph_settings["show_grid"] = st.checkbox(
                    "Show Grid", value=True
                )
                st.session_state.graph_settings["bar_gap"] = st.slider(
                    "Bar Gap", 0.0, 0.5, 0.15
                )
            with c3:
                st.session_state.graph_settings["font_size"] = st.slider(
                    "Font Size", 10, 24, 14
                )
                st.session_state.graph_settings["figure_height"] = st.slider(
                    "Height", 400, 1200, 600
                )
                st.session_state.graph_settings["figure_width"] = st.slider(
                    "Width", 600, 1200, 900
                )

        st.divider()

        # --- Data preparation ---
        data = st.session_state.processed_data

        # Handle both dict-of-genes or single DataFrame
        if isinstance(data, dict):
            all_genes = list(data.keys())
        else:
            if "Target" in data.columns:
                all_genes = sorted(data["Target"].dropna().unique().tolist())
            else:
                st.error("No 'Target' column found in processed data.")
                st.stop()

        selected_genes = st.multiselect(
            "🧬 Select Genes to Visualize",
            all_genes,
            default=all_genes[: min(3, len(all_genes))],
        )

        st.divider()

        # --- Draw selected gene graphs ---
        for gene in selected_genes:
            st.markdown(f"### **{gene}**")

            # get data per gene
            if isinstance(data, dict):
                df = data[gene]
            else:
                df = data[data["Target"] == gene].copy()

            if df.empty:
                st.warning(f"No data available for {gene}")
                continue

            # choose conditions to display
            conds = df["Condition"].unique().tolist()
            visible = [
                c
                for c in conds
                if st.checkbox(f"Show {c}", value=True, key=f"chk_{gene}_{c}")
            ]
            df_filtered = df[df["Condition"].isin(visible)]

            # choose color
            if gene not in st.session_state.graph_settings["bar_colors"]:
                st.session_state.graph_settings["bar_colors"][gene] = px.colors.qualitative.Plotly[
                    len(st.session_state.graph_settings["bar_colors"]) % 10
                ]
            color = st.color_picker(
                "Bar Color",
                st.session_state.graph_settings["bar_colors"][gene],
                key=f"color_{gene}",
            )
            st.session_state.graph_settings["bar_colors"][gene] = color

            # create figure
            fig = GraphGenerator.create_gene_graph(
                df_filtered,
                gene,
                st.session_state.graph_settings,
                sample_order=st.session_state.get("sample_order", None),
            )

            st.plotly_chart(fig, use_container_width=True, key=f"fig_{gene}")
            st.divider()

# ==================== TAB 4: GRAPHS ====================
with tab4:
    st.header("Step 4: Individual Gene Graphs")
    
    if st.session_state.processed_data:
        # Initialize graph settings (extend with new options)
        if 'graph_settings' not in st.session_state:
            st.session_state.graph_settings = {
                'title_size': 20, 'font_size': 14, 'sig_font_size': 16,
                'figure_width': 1000, 'figure_height': 600,
                'color_scheme': 'plotly_white', 'show_error': True,
                'show_significance': True, 'show_grid': True,
                'xlabel': 'Condition', 'ylabel': 'Fold Change (Relative to Control)',
                'bar_colors': {},
                # NEW options
                'orientation': 'v',  # 'v' or 'h'
                'error_multiplier': 1.96,
                'bar_opacity': 0.95,
                'bar_gap': 0.15,
                'marker_line_width': 1,
                'show_legend': False,
                'y_log_scale': False,
                'y_min': None,
                'y_max': None
            }

        # Global settings
        st.subheader("⚙️ Global Graph Settings")

        # layout widgets in groups for clarity
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.session_state.graph_settings['font_size'] = st.slider("Font Size", 8, 28, st.session_state.graph_settings['font_size'])
            st.session_state.graph_settings['title_size'] = st.slider("Title Size", 10, 36, st.session_state.graph_settings['title_size'])
            st.session_state.graph_settings['sig_font_size'] = st.slider("Significance Font Size", 8, 28, st.session_state.graph_settings['sig_font_size'])
        with col2:
            st.session_state.graph_settings['figure_width'] = st.slider("Width (px)", 600, 2000, st.session_state.graph_settings['figure_width'])
            st.session_state.graph_settings['figure_height'] = st.slider("Height (px)", 300, 1400, st.session_state.graph_settings['figure_height'])
            st.session_state.graph_settings['orientation'] = st.selectbox("Bar Orientation", ['v', 'h'], index=0 if st.session_state.graph_settings['orientation']=='v' else 1)
        with col3:
            st.session_state.graph_settings['color_scheme'] = st.selectbox(
                "Theme", ['plotly_white', 'plotly', 'plotly_dark', 'seaborn', 'simple_white', 'presentation'],
                index=0 if st.session_state.graph_settings.get('color_scheme','plotly_white')=='plotly_white' else 0
            )
            st.session_state.graph_settings['show_legend'] = st.checkbox("Show Legend", st.session_state.graph_settings.get('show_legend', False))
            st.session_state.graph_settings['show_grid'] = st.checkbox("Grid", st.session_state.graph_settings.get('show_grid', True))
        with col4:
            st.session_state.graph_settings['show_error'] = st.checkbox("Error Bars", st.session_state.graph_settings.get('show_error', True))
            st.session_state.graph_settings['show_significance'] = st.checkbox("Significance", st.session_state.graph_settings.get('show_significance', True))
            st.session_state.graph_settings['error_multiplier'] = st.slider("Error multiplier (σ → CI)", 0.5, 3.0, float(st.session_state.graph_settings.get('error_multiplier',1.96)), step=0.01)

        # more detailed controls
        col5, col6, col7 = st.columns([1.2, 1.2, 1])
        with col5:
            st.session_state.graph_settings['bar_opacity'] = st.slider("Bar Opacity", 0.1, 1.0, st.session_state.graph_settings.get('bar_opacity',0.95))
            st.session_state.graph_settings['bar_gap'] = st.slider("Bar Gap (group gap)", 0.0, 0.5, st.session_state.graph_settings.get('bar_gap',0.15))
        with col6:
            st.session_state.graph_settings['marker_line_width'] = st.slider("Bar Edge Width", 0, 6, st.session_state.graph_settings.get('marker_line_width',1))
            st.session_state.graph_settings['y_log_scale'] = st.checkbox("Y-axis log scale", st.session_state.graph_settings.get('y_log_scale', False))
        with col7:
            y_range = st.slider("Y-axis range (min,max) - set both to adjust, leave unchanged for auto", 0.0, 100.0, (0.0, 10.0))
            # if user hasn't changed from defaults, keep None to use auto-scaling
            if y_range != (0.0, 10.0):
                st.session_state.graph_settings['y_min'] = float(y_range[0])
                st.session_state.graph_settings['y_max'] = float(y_range[1])
            else:
                st.session_state.graph_settings['y_min'] = None
                st.session_state.graph_settings['y_max'] = None

        # --- Minimal preset controls removed (user requested elimination of presets) ---
        col1, col2 = st.columns([1,1])
        with col1:
            if st.button("🔄 Reset Graph Settings"):
                if 'graph_settings' in st.session_state:
                    del st.session_state['graph_settings']
                st.rerun()
        with col2:
            st.markdown("Presets removed — use sliders & controls to customize.")

        st.markdown("---")

        # Generate graphs for each gene
        st.subheader("📊 Gene-Specific Graphs")
        
        efficacy_config = EFFICACY_CONFIG.get(st.session_state.selected_efficacy, {})
        
        # Create graphs
        if 'graphs' not in st.session_state:
            st.session_state.graphs = {}
        
        for gene in st.session_state.processed_data.keys():
            st.markdown(f"### 🧬 {gene}")
            
            col_left, col_right = st.columns([1.2, 4])
            with col_left:
                # gene-level color
                st.session_state.graph_settings['bar_colors'].setdefault(gene, px.colors.qualitative.Plotly[len(st.session_state.graph_settings['bar_colors']) % len(px.colors.qualitative.Plotly)])
                st.session_state.graph_settings['bar_colors'][gene] = st.color_picker(f"{gene} color", st.session_state.graph_settings['bar_colors'][gene], key=f"color_{gene}")

                # per-sample per-gene override store
                if 'sample_gene_overrides' not in st.session_state:
                    st.session_state['sample_gene_overrides'] = {}  # key = (gene, sample) -> dict{'include':bool, 'color':hex or None}
                st.markdown("**Per-sample overrides**")
                # list conditions for this gene (from processed_data)
                conds = st.session_state.processed_data[gene]['Condition'].tolist()
                for cond in conds:
                    key = (gene, cond)
                    # default include True unless overridden
                    cur = st.session_state.sample_gene_overrides.get(str(key), {'include': True, 'color': None})
                    inc = st.checkbox(f"Include {cond}", value=cur.get('include', True), key=f"incl_{gene}_{cond}")
                    col = st.text_input(f"Color for {cond} (hex or blank)", value=cur.get('color') or '', key=f"col_{gene}_{cond}", label_visibility="collapsed")
                    st.session_state.sample_gene_overrides[str(key)] = {'include': bool(inc), 'color': col.strip() or None}

            with col_right:
                # pass overrides into GraphGenerator via the settings object (non-destructive)
                local_settings = dict(st.session_state.graph_settings)
                local_settings['bar_colors'] = local_settings.get('bar_colors', {}).copy()
                # apply per-sample color overrides (for this gene)
                for cond in conds:
                    ov = st.session_state.sample_gene_overrides.get(str((gene, cond)), {})
                    if ov.get('color'):
                        local_settings['bar_colors'][f"{gene}__{cond}"] = ov['color']

        @staticmethod
        def create_gene_graph(data: pd.DataFrame, gene: str, settings: dict, efficacy_config: dict = None, sample_order: List[str] = None, per_sample_overrides: dict = None) -> go.Figure:
            """Create individual graph for each gene with support for sample ordering and per-sample overrides"""
            gene_data = data[data['Target'] == gene].copy()

            # Map group order and then apply a user-specified sample order if provided
            group_order = ['Baseline', 'Negative Control', 'Positive Control', 'Treatment']
            gene_data['group_sort'] = gene_data['Group'].apply(lambda x: group_order.index(x) if x in group_order else 999)

            # Apply sample_order if provided to set custom Condition order (filter to included ones)
            if sample_order:
                # Keep only conditions in sample_order (but if gene has others, append them)
                ordered = [c for c in sample_order if c in gene_data['Condition'].unique()]
                other = [c for c in sorted(gene_data['Condition'].unique()) if c not in ordered]
                final_order = ordered + other
                gene_data['Condition'] = pd.Categorical(gene_data['Condition'], categories=final_order, ordered=True)
                gene_data = gene_data.sort_values(['group_sort','Condition'])
            else:
                gene_data = gene_data.sort_values(['group_sort', 'Condition'])

            # Apply per-sample-per-gene include overrides: drop excluded ones
            if per_sample_overrides:
                # keys are stored as str((gene, condition))
                keep_mask = []
                for _, row in gene_data.iterrows():
                    key = str((gene, row['Condition']))
                    ov = per_sample_overrides.get(key, {})
                    # default include True
                    keep_mask.append(bool(ov.get('include', True)))
                gene_data = gene_data[np.array(keep_mask)]

            fig = go.Figure()

            # Color mapping by group (and fallbacks)
            default_color = settings.get('bar_colors', {}).get(gene, '#95E1D3')
            color_map = {
                'Baseline': '#808080',
                'Negative Control': '#FF6B6B',
                'Positive Control': '#4ECDC4',
                'Treatment': default_color
            }

            # per-bar colors: respect per-sample overrides if present
            bar_colors = []
            for _, row in gene_data.iterrows():
                key = str((gene, row['Condition']))
                ov = (per_sample_overrides or {}).get(key, {})
                if ov and ov.get('color'):
                    bar_colors.append(ov.get('color'))
                else:
                    bar_colors.append(color_map.get(row['Group'], default_color))

            # error bars (apply multiplier option)
            err_array = (gene_data['SEM'] * settings.get('error_multiplier', 1.96)).fillna(0).values

            fig.add_trace(go.Bar(
                x=gene_data['Condition'],
                y=gene_data['Relative_Expression'],
                error_y=dict(type='data', array=err_array, visible=settings.get('show_error', True)),
                text=gene_data['significance'] if settings.get('show_significance', True) else None,
                textposition='outside',
                textfont=dict(size=settings.get('sig_font_size', 16)),
                marker=dict(color=bar_colors, line=dict(width=settings.get('marker_line_width', 1), color='black')),
                showlegend=settings.get('show_legend', False),
                opacity=settings.get('bar_opacity', 0.95)
            ))

            # Add expected direction indicator
            if efficacy_config and 'expected_direction' in efficacy_config:
                direction = efficacy_config.get('expected_direction', {}).get(gene, '')
                direction_text = '↑ Expected increase' if direction == 'up' else '↓ Expected decrease' if direction == 'down' else ''
                if direction_text:
                    fig.add_annotation(text=direction_text, xref='paper', yref='paper', x=0.02, y=0.98, showarrow=False, font=dict(size=12, color='red'), align='left')

            # Layout
            layout = dict(
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
            if settings.get('y_log_scale'):
                layout['yaxis']['type'] = 'log'
            if settings.get('y_min') is not None or settings.get('y_max') is not None:
                ymin = settings.get('y_min')
                ymax = settings.get('y_max')
                # Plotly expects range as [min,max]; None values are ignored
                layout['yaxis']['range'] = [ymin if ymin is not None else None, ymax if ymax is not None else None]

            layout['bargap'] = settings.get('bar_gap', 0.15)

            fig.update_layout(**layout)
            return fig
        
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
