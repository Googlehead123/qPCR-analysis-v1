import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from scipy import stats
import io
import json
from datetime import datetime
import re

# Page config
st.set_page_config(page_title="qPCR Analysis Suite", layout="wide", initial_sidebar_state="expanded")

# Initialize session state
if 'data' not in st.session_state:
    st.session_state.data = None
if 'processed_data' not in st.session_state:
    st.session_state.processed_data = None
if 'sample_mapping' not in st.session_state:
    st.session_state.sample_mapping = {}
if 'analysis_templates' not in st.session_state:
    st.session_state.analysis_templates = {}
if 'current_graph' not in st.session_state:
    st.session_state.current_graph = None
if 'excluded_wells' not in st.session_state:
    st.session_state.excluded_wells = set()

# Efficacy database
EFFICACY_DB = {
    '항노화': {'genes': ['MMP-1'], 'cell': 'HS68', 'controls': 1},
    '탄력': {'genes': ['COL1A1', 'ELN', 'FBN-1'], 'cell': 'HS68', 'controls': 1},
    '재생': {'genes': ['Migration'], 'cell': 'HS68', 'controls': 1},
    '보습': {'genes': ['HAS3', 'AQP3'], 'cell': 'HaCaT', 'controls': 1},
    '장벽': {'genes': ['FLG', 'CLDN', 'IVL'], 'cell': 'HaCaT', 'controls': 1},
    '진정': {'genes': ['IL-1β', 'IL-6', 'TNFα'], 'cell': 'HaCaT', 'controls': 2},
    '가려움': {'genes': ['TSLP'], 'cell': 'HaCaT', 'controls': 1},
    '표피증식': {'genes': ['PCNA', 'KI67'], 'cell': 'HaCaT', 'controls': 1},
    '멜라닌억제': {'genes': ['MITF', 'TYR'], 'cell': 'B16F10', 'controls': 1},
    '지질억제': {'genes': ['SREBPa', 'SREBPc'], 'cell': 'SZ95', 'controls': 1},
}

class QPCRParser:
    """Smart parser for multiple qPCR data formats"""
    
    @staticmethod
    def detect_format(df):
        """Auto-detect data format"""
        if 'Well Position' in df.columns:
            return 'format1'
        elif 'Sample Name' in df.columns and 'Target Name' in df.columns:
            return 'format2'
        return 'unknown'
    
    @staticmethod
    def parse_format1(df):
        """Parse QuantStudio format"""
        # Skip metadata rows
        data_start = df[df.iloc[:, 0] == 'Well'].index[0]
        df = df.iloc[data_start:].reset_index(drop=True)
        df.columns = df.iloc[0]
        df = df.iloc[1:].reset_index(drop=True)
        
        # Extract relevant columns
        parsed = pd.DataFrame({
            'Well': df['Well Position'],
            'Sample': df['Sample Name'],
            'Target': df['Target Name'],
            'CT': pd.to_numeric(df['CT'], errors='coerce'),
            'CT_Mean': pd.to_numeric(df['Ct Mean'], errors='coerce'),
            'CT_SD': pd.to_numeric(df['Ct SD'], errors='coerce'),
            'RQ': pd.to_numeric(df['RQ'], errors='coerce'),
            'Delta_CT': pd.to_numeric(df['Delta Ct'], errors='coerce'),
            'Delta_CT_Mean': pd.to_numeric(df['Delta Ct Mean'], errors='coerce'),
        })
        return parsed.dropna(subset=['CT'])
    
    @staticmethod
    def parse_format2(df):
        """Parse StepOnePlus format"""
        # Skip metadata rows
        data_start = df[df.iloc[:, 0] == 'Well'].index[0]
        df = df.iloc[data_start:].reset_index(drop=True)
        df.columns = df.iloc[0]
        df = df.iloc[1:].reset_index(drop=True)
        
        parsed = pd.DataFrame({
            'Well': df['Well'],
            'Sample': df['Sample Name'],
            'Target': df['Target Name'],
            'CT': pd.to_numeric(df['Cт'], errors='coerce'),
            'CT_Mean': pd.to_numeric(df['Cт Mean'], errors='coerce'),
            'CT_SD': pd.to_numeric(df['Cт SD'], errors='coerce'),
            'RQ': pd.to_numeric(df['RQ'], errors='coerce'),
            'Delta_CT': pd.to_numeric(df['ΔCт'], errors='coerce'),
            'Delta_CT_Mean': pd.to_numeric(df['ΔCт Mean'], errors='coerce'),
        })
        return parsed.dropna(subset=['CT'])
    
    @staticmethod
    def parse(file):
        """Main parsing function"""
        df = pd.read_csv(file, encoding='utf-8', low_memory=False)
        format_type = QPCRParser.detect_format(df)
        
        if format_type == 'format1':
            return QPCRParser.parse_format1(df)
        elif format_type == 'format2':
            return QPCRParser.parse_format2(df)
        else:
            st.error("Unknown format. Please check your file.")
            return None

class StatisticalAnalyzer:
    """Statistical analysis engine"""
    
    @staticmethod
    def calculate_ddct(data, reference_sample, reference_gene, excluded_wells):
        """Calculate ΔΔCt and relative expression"""
        results = []
        
        # Filter out excluded wells
        data = data[~data['Well'].isin(excluded_wells)].copy()
        
        # Get reference values
        ref_data = data[(data['Sample'] == reference_sample) & (data['Target'] == reference_gene)]
        if len(ref_data) == 0:
            st.error(f"Reference sample '{reference_sample}' or gene '{reference_gene}' not found")
            return None
        
        ref_ct_mean = ref_data['CT'].mean()
        
        # Calculate for each sample-target combination
        for sample in data['Sample'].unique():
            for target in data['Target'].unique():
                subset = data[(data['Sample'] == sample) & (data['Target'] == target)]
                if len(subset) == 0:
                    continue
                
                # Get housekeeping gene CT for this sample
                hk_data = data[(data['Sample'] == sample) & (data['Target'] == reference_gene)]
                if len(hk_data) == 0:
                    continue
                
                hk_ct_mean = hk_data['CT'].mean()
                target_ct_mean = subset['CT'].mean()
                target_ct_sd = subset['CT'].std()
                
                # ΔCt = Ct(target) - Ct(housekeeping)
                delta_ct = target_ct_mean - hk_ct_mean
                
                # ΔΔCt = ΔCt(sample) - ΔCt(reference)
                ref_hk_ct = ref_data['CT'].mean()
                ref_target_data = data[(data['Sample'] == reference_sample) & (data['Target'] == target)]
                if len(ref_target_data) > 0:
                    ref_target_ct = ref_target_data['CT'].mean()
                    ref_delta_ct = ref_target_ct - ref_hk_ct
                else:
                    ref_delta_ct = 0
                
                ddct = delta_ct - ref_delta_ct
                
                # Relative expression = 2^(-ΔΔCt)
                rel_expr = 2 ** (-ddct)
                
                results.append({
                    'Sample': sample,
                    'Target': target,
                    'CT_Mean': target_ct_mean,
                    'CT_SD': target_ct_sd,
                    'n_replicates': len(subset),
                    'Delta_CT': delta_ct,
                    'Delta_Delta_CT': ddct,
                    'Relative_Expression': rel_expr,
                    'SEM': target_ct_sd / np.sqrt(len(subset)) if len(subset) > 1 else 0
                })
        
        return pd.DataFrame(results)
    
    @staticmethod
    def perform_statistics(processed_data, control_sample):
        """Perform t-tests comparing each sample to control"""
        results = processed_data.copy()
        results['p_value'] = np.nan
        results['significance'] = ''
        
        for target in results['Target'].unique():
            control_expr = results[(results['Sample'] == control_sample) & 
                                  (results['Target'] == target)]['Relative_Expression'].values
            
            for sample in results['Sample'].unique():
                if sample == control_sample:
                    continue
                
                sample_expr = results[(results['Sample'] == sample) & 
                                     (results['Target'] == target)]['Relative_Expression'].values
                
                if len(control_expr) > 0 and len(sample_expr) > 0:
                    # Perform t-test (assuming we have raw replicates, otherwise use means)
                    t_stat, p_val = stats.ttest_ind(control_expr, sample_expr)
                    
                    idx = (results['Sample'] == sample) & (results['Target'] == target)
                    results.loc[idx, 'p_value'] = p_val
                    
                    # Significance stars
                    if p_val < 0.001:
                        results.loc[idx, 'significance'] = '***'
                    elif p_val < 0.01:
                        results.loc[idx, 'significance'] = '**'
                    elif p_val < 0.05:
                        results.loc[idx, 'significance'] = '*'
        
        return results

class GraphGenerator:
    """Publication-ready graph generator"""
    
    @staticmethod
    def create_bar_chart(data, title, ylabel, show_error_bars=True):
        """Create grouped bar chart with error bars"""
        fig = go.Figure()
        
        samples = data['Sample'].unique()
        targets = data['Target'].unique()
        
        for target in targets:
            target_data = data[data['Target'] == target]
            
            fig.add_trace(go.Bar(
                name=target,
                x=target_data['Sample'],
                y=target_data['Relative_Expression'],
                error_y=dict(
                    type='data',
                    array=target_data['SEM'] * 1.96 if show_error_bars else None,  # 95% CI
                    visible=show_error_bars
                ),
                text=target_data['significance'],
                textposition='outside',
                textfont=dict(size=16)
            ))
        
        fig.update_layout(
            title=dict(text=title, font=dict(size=20)),
            xaxis_title="Sample",
            yaxis_title=ylabel,
            barmode='group',
            template='plotly_white',
            font=dict(size=14),
            height=600,
            showlegend=True,
            legend=dict(font=dict(size=14))
        )
        
        return fig
    
    @staticmethod
    def create_heatmap(data):
        """Create heatmap for multi-gene analysis"""
        pivot = data.pivot(index='Target', columns='Sample', values='Relative_Expression')
        
        fig = go.Figure(data=go.Heatmap(
            z=pivot.values,
            x=pivot.columns,
            y=pivot.index,
            colorscale='RdYlGn',
            text=np.round(pivot.values, 2),
            texttemplate='%{text}',
            textfont={"size": 12},
            colorbar=dict(title="Fold Change")
        ))
        
        fig.update_layout(
            title="Gene Expression Heatmap",
            xaxis_title="Sample",
            yaxis_title="Gene",
            height=400 + len(pivot.index) * 30,
            font=dict(size=14)
        )
        
        return fig

def export_to_excel(raw_data, processed_data, analysis_params):
    """Export comprehensive Excel file with all calculations"""
    output = io.BytesIO()
    
    with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
        # Sheet 1: Analysis parameters
        params_df = pd.DataFrame([analysis_params])
        params_df.to_excel(writer, sheet_name='Analysis_Parameters', index=False)
        
        # Sheet 2: Raw data
        raw_data.to_excel(writer, sheet_name='Raw_Data', index=False)
        
        # Sheet 3: Processed data with calculations
        processed_data.to_excel(writer, sheet_name='Calculations', index=False)
        
        # Sheet 4: Summary statistics
        summary = processed_data.groupby('Target').agg({
            'Relative_Expression': ['mean', 'std', 'sem'],
            'p_value': 'min'
        }).round(4)
        summary.to_excel(writer, sheet_name='Summary')
    
    return output.getvalue()

# ==================== STREAMLIT UI ====================

st.title("🧬 qPCR Analysis Automation Suite")
st.markdown("**Automated analysis from raw data to publication-ready graphs**")

# Sidebar - Chat Assistant
with st.sidebar:
    st.header("💬 Quick Commands")
    st.markdown("""
    **Tips for faster analysis:**
    - Upload multiple files at once
    - Use templates for repeated assays
    - Flag outliers instead of deleting
    - Save templates for team reuse
    """)
    
    # Template management
    st.subheader("📋 Analysis Templates")
    template_name = st.text_input("Save current analysis as:")
    if st.button("💾 Save Template") and st.session_state.sample_mapping:
        st.session_state.analysis_templates[template_name] = {
            'mapping': st.session_state.sample_mapping.copy(),
            'timestamp': datetime.now().isoformat()
        }
        st.success(f"Template '{template_name}' saved!")
    
    if st.session_state.analysis_templates:
        load_template = st.selectbox("Load template:", [""] + list(st.session_state.analysis_templates.keys()))
        if load_template:
            st.session_state.sample_mapping = st.session_state.analysis_templates[load_template]['mapping']
            st.info(f"Template '{load_template}' loaded!")

# Main content tabs
tab1, tab2, tab3, tab4, tab5 = st.tabs(["📁 Upload", "🗺️ Sample Mapping", "🔬 Analysis", "📊 Visualization", "📤 Export"])

# TAB 1: File Upload
with tab1:
    st.header("Step 1: Upload qPCR Data")
    
    uploaded_files = st.file_uploader(
        "Upload CSV file(s) from QuantStudio or StepOnePlus",
        type=['csv'],
        accept_multiple_files=True
    )
    
    if uploaded_files:
        all_data = []
        for file in uploaded_files:
            parsed = QPCRParser.parse(file)
            if parsed is not None:
                parsed['Source_File'] = file.name
                all_data.append(parsed)
                st.success(f"✅ Parsed {file.name}: {len(parsed)} wells")
        
        if all_data:
            st.session_state.data = pd.concat(all_data, ignore_index=True)
            
            # Display preview
            st.subheader("Data Preview")
            col1, col2, col3 = st.columns(3)
            col1.metric("Total Wells", len(st.session_state.data))
            col2.metric("Samples", st.session_state.data['Sample'].nunique())
            col3.metric("Genes", st.session_state.data['Target'].nunique())
            
            # Quality metrics
            st.subheader("Quality Metrics")
            qc_data = st.session_state.data.groupby('Target').agg({
                'CT': ['mean', 'std', 'min', 'max']
            }).round(2)
            st.dataframe(qc_data)
            
            # Flag outliers
            st.subheader("🚩 Outlier Detection")
            for target in st.session_state.data['Target'].unique():
                target_data = st.session_state.data[st.session_state.data['Target'] == target]
                q1 = target_data['CT'].quantile(0.25)
                q3 = target_data['CT'].quantile(0.75)
                iqr = q3 - q1
                outliers = target_data[(target_data['CT'] < q1 - 1.5*iqr) | (target_data['CT'] > q3 + 1.5*iqr)]
                
                if len(outliers) > 0:
                    st.warning(f"**{target}**: {len(outliers)} potential outliers detected")
                    with st.expander(f"View {target} outliers"):
                        st.dataframe(outliers[['Well', 'Sample', 'CT']])
                        if st.checkbox(f"Exclude these wells from {target} analysis", key=f"exclude_{target}"):
                            st.session_state.excluded_wells.update(outliers['Well'].tolist())
            
            st.dataframe(st.session_state.data.head(20))

# TAB 2: Sample Mapping
with tab2:
    st.header("Step 2: Map Sample Names to Conditions")
    
    if st.session_state.data is not None:
        st.markdown("**Assign biological meaning to sample labels**")
        
        # Auto-detect efficacy type
        detected_genes = st.session_state.data['Target'].unique()
        suggested_efficacy = None
        for efficacy, info in EFFICACY_DB.items():
            if any(gene in detected_genes for gene in info['genes']):
                suggested_efficacy = efficacy
                break
        
        col1, col2 = st.columns(2)
        with col1:
            efficacy_type = st.selectbox(
                "Efficacy Test Type",
                options=list(EFFICACY_DB.keys()),
                index=list(EFFICACY_DB.keys()).index(suggested_efficacy) if suggested_efficacy else 0
            )
        
        with col2:
            expected_genes = EFFICACY_DB[efficacy_type]['genes']
            st.info(f"Expected genes: {', '.join(expected_genes)}")
        
        # Sample mapping interface
        st.subheader("Sample Mapping")
        samples = sorted(st.session_state.data['Sample'].unique(), key=lambda x: (not str(x).isdigit(), str(x)))
        
        mapping_data = []
        for sample in samples:
            if sample not in st.session_state.sample_mapping:
                st.session_state.sample_mapping[sample] = {
                    'condition': sample,
                    'group': 'Treatment',
                    'concentration': ''
                }
            
            col1, col2, col3, col4 = st.columns([2, 3, 2, 2])
            with col1:
                st.text(sample)
            with col2:
                condition = st.text_input(
                    "Condition Name",
                    value=st.session_state.sample_mapping[sample]['condition'],
                    key=f"cond_{sample}"
                )
                st.session_state.sample_mapping[sample]['condition'] = condition
            with col3:
                group = st.selectbox(
                    "Group",
                    ['Control', 'Positive Control', 'Treatment', 'Induced Control'],
                    index=['Control', 'Positive Control', 'Treatment', 'Induced Control'].index(
                        st.session_state.sample_mapping[sample]['group']
                    ),
                    key=f"group_{sample}"
                )
                st.session_state.sample_mapping[sample]['group'] = group
            with col4:
                conc = st.text_input(
                    "Concentration",
                    value=st.session_state.sample_mapping[sample]['concentration'],
                    key=f"conc_{sample}"
                )
                st.session_state.sample_mapping[sample]['concentration'] = conc
        
        # Show mapping summary
        st.subheader("Mapping Summary")
        mapping_df = pd.DataFrame([
            {'Original': k, **v} for k, v in st.session_state.sample_mapping.items()
        ])
        st.dataframe(mapping_df)
    else:
        st.warning("⚠️ Please upload data first")

# TAB 3: Analysis
with tab3:
    st.header("Step 3: Statistical Analysis")
    
    if st.session_state.data is not None and st.session_state.sample_mapping:
        # Select reference sample and gene
        col1, col2 = st.columns(2)
        
        with col1:
            control_samples = [k for k, v in st.session_state.sample_mapping.items() 
                             if v['group'] == 'Control']
            reference_sample = st.selectbox("Reference Sample (Negative Control)", control_samples)
        
        with col2:
            housekeeping_genes = [g for g in st.session_state.data['Target'].unique() 
                                 if g.upper() in ['ACTIN', 'B-ACTIN', 'GAPDH', 'ACTB']]
            reference_gene = st.selectbox("Housekeeping Gene", housekeeping_genes)
        
        if st.button("🔬 Run Analysis", type="primary"):
            with st.spinner("Calculating ΔΔCt and statistics..."):
                # Calculate ΔΔCt
                processed = StatisticalAnalyzer.calculate_ddct(
                    st.session_state.data,
                    reference_sample,
                    reference_gene,
                    st.session_state.excluded_wells
                )
                
                if processed is not None:
                    # Perform statistics
                    processed = StatisticalAnalyzer.perform_statistics(processed, reference_sample)
                    
                    # Apply sample mapping
                    processed['Condition'] = processed['Sample'].map(
                        lambda x: st.session_state.sample_mapping.get(x, {}).get('condition', x)
                    )
                    processed['Group'] = processed['Sample'].map(
                        lambda x: st.session_state.sample_mapping.get(x, {}).get('group', 'Unknown')
                    )
                    
                    st.session_state.processed_data = processed
                    st.success("✅ Analysis complete!")
        
        # Display results
        if st.session_state.processed_data is not None:
            st.subheader("Analysis Results")
            
            # Summary statistics
            col1, col2, col3 = st.columns(3)
            sig_count = (st.session_state.processed_data['p_value'] < 0.05).sum()
            col1.metric("Significant Results", f"{sig_count} / {len(st.session_state.processed_data)}")
            col2.metric("Mean Fold Change", 
                       f"{st.session_state.processed_data['Relative_Expression'].mean():.2f}x")
            col3.metric("Max Fold Change", 
                       f"{st.session_state.processed_data['Relative_Expression'].max():.2f}x")
            
            # Data table with filtering
            st.subheader("Detailed Results")
            show_filter = st.checkbox("Show only significant results (p < 0.05)")
            
            display_data = st.session_state.processed_data.copy()
            if show_filter:
                display_data = display_data[display_data['p_value'] < 0.05]
            
            st.dataframe(
                display_data.style.background_gradient(subset=['Relative_Expression'], cmap='RdYlGn')
                .format({'Relative_Expression': '{:.3f}', 'p_value': '{:.4f}', 'CT_Mean': '{:.2f}'})
            )
    else:
        st.warning("⚠️ Complete previous steps first")

# TAB 4: Visualization
with tab4:
    st.header("Step 4: Generate Publication-Ready Graphs")
    
    if st.session_state.processed_data is not None:
        # Graph customization
        col1, col2 = st.columns(2)
        with col1:
            graph_type = st.selectbox("Graph Type", ['Bar Chart', 'Heatmap', 'Both'])
            title = st.text_input("Graph Title", "Relative Gene Expression")
        
        with col2:
            ylabel = st.text_input("Y-axis Label", "Fold Change (Relative to Control)")
            show_error = st.checkbox("Show Error Bars", value=True)
        
        # Generate graph
        if graph_type in ['Bar Chart', 'Both']:
            st.subheader("Bar Chart")
            fig_bar = GraphGenerator.create_bar_chart(
                st.session_state.processed_data,
                title,
                ylabel,
                show_error
            )
            st.plotly_chart(fig_bar, use_container_width=True)
            st.session_state.current_graph = fig_bar
        
        if graph_type in ['Heatmap', 'Both']:
            st.subheader("Heatmap")
            fig_heat = GraphGenerator.create_heatmap(st.session_state.processed_data)
            st.plotly_chart(fig_heat, use_container_width=True)
        
        # Interactive editing
        st.subheader("🎨 Interactive Graph Editor")
        with st.expander("Customize appearance"):
            col1, col2, col3 = st.columns(3)
            with col1:
                font_size = st.slider("Font Size", 10, 24, 14)
            with col2:
                bar_width = st.slider("Bar Width", 0.3, 1.0, 0.8)
            with col3:
                color_scheme = st.selectbox("Color Scheme", 
                    ['plotly', 'ggplot2', 'seaborn', 'simple_white'])
            
            if st.button("Apply Changes"):
                fig_bar.update_layout(font=dict(size=font_size), template=color_scheme)
                st.plotly_chart(fig_bar, use_container_width=True)
    else:
        st.warning("⚠️ Run analysis first")

# TAB 5: Export
with tab5:
    st.header("Step 5: Export Results")
    
    if st.session_state.processed_data is not None:
        st.subheader("Download Options")
        
        col1, col2 = st.columns(2)
        
        with col1:
            # Excel export
            st.markdown("### 📊 Excel (Full Analysis)")
            analysis_params = {
                'Date': datetime.now().strftime("%Y-%m-%d %H:%M"),
                'Reference_Sample': reference_sample if 'reference_sample' in locals() else 'N/A',
                'Housekeeping_Gene': reference_gene if 'reference_gene' in locals() else 'N/A',
                'Excluded_Wells': len(st.session_state.excluded_wells)
            }
            
            excel_data = export_to_excel(
                st.session_state.data,
                st.session_state.processed_data,
                analysis_params
            )
            
            st.download_button(
                label="📥 Download Excel Report",
                data=excel_data,
                file_name=f"qPCR_analysis_{datetime.now().strftime('%Y%m%d_%H%M')}.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )
        
        with col2:
            # Graph export
            st.markdown("### 📈 Graph (PNG/SVG)")
            if st.session_state.current_graph:
                img_format = st.radio("Format", ['PNG', 'SVG'])
                
                img_bytes = st.session_state.current_graph.to_image(
                    format=img_format.lower(),
                    width=1200,
                    height=800,
                    scale=2
                )
                
                st.download_button(
                    label=f"📥 Download {img_format}",
                    data=img_bytes,
                    file_name=f"qPCR_graph_{datetime.now().strftime('%Y%m%d_%H%M')}.{img_format.lower()}",
                    mime=f"image/{img_format.lower()}"
                )
        
        # JSON export for reproducibility
        st.subheader("🔄 Reproducibility Package")
        repro_data = {
            'analysis_date': datetime.now().isoformat(),
            'sample_mapping': st.session_state.sample_mapping,
            'excluded_wells': list(st.session_state.excluded_wells),
            'parameters': analysis_params if 'analysis_params' in locals() else {}
        }
        
        st.download_button(
            label="📥 Download Analysis Config (JSON)",
            data=json.dumps(repro_data, indent=2),
            file_name=f"qPCR_config_{datetime.now().strftime('%Y%m%d_%H%M')}.json",
            mime="application/json"
        )
        
        st.success("✅ All export options ready!")
    else:
        st.warning("⚠️ Complete analysis first")

# Footer
st.markdown("---")
st.markdown("""
<div style='text-align: center; color: #666;'>
    <p>🧬 qPCR Analysis Suite v1.0 | Built for precision cosmetic science research</p>
    <p>Save templates • Flag outliers • Generate publication-ready graphs</p>
</div>
""", unsafe_allow_html=True)