import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from scipy import stats
from difflib import SequenceMatcher
import io
import json
from datetime import datetime
from typing import Dict, List, Optional
import re

# ==================== CONFIG ====================
st.set_page_config(page_title="qPCR Analysis Ultimate", layout="wide")

# Session state
defaults = {
    'data': None, 'processed_data': {}, 'sample_mapping': {}, 'analysis_templates': {},
    'graphs': {}, 'excluded_wells': set(), 'excluded_samples': set(), 'selected_efficacy': None,
    'hk_gene': None, 'included_samples_analysis': {}, 'graph_settings': {}, 'auto_merged': {}
}
for k, v in defaults.items():
    if k not in st.session_state:
        st.session_state[k] = v

# Efficacy configs
EFFICACY_CONFIG = {
    '탄력': {'genes': ['COL1A1', 'ELN', 'FBN-1', 'FBN1'], 'cell': 'HS68'},
    '항노화': {'genes': ['COL1A1', 'COL1', 'MMP-1', 'MMP1'], 'cell': 'HS68',
              'expected': {'COL1A1': '↑', 'MMP-1': '↓'}},
    '보습': {'genes': ['AQP3', 'HAS3'], 'cell': 'HaCaT'},
    '장벽': {'genes': ['FLG', 'CLDN', 'IVL'], 'cell': 'HaCaT'},
    '표피증식': {'genes': ['KI67', 'PCNA'], 'cell': 'HaCaT'},
    '멜라닌억제': {'genes': ['MITF', 'TYR'], 'cell': 'B16F10', 'expected': {'MITF': '↓', 'TYR': '↓'}},
    '진정': {'genes': ['IL1B', 'IL-1β', 'IL6', 'TNFA', 'TNFα'], 'cell': 'HaCaT',
           'expected': {'IL1B': '↓', 'IL6': '↓', 'TNFA': '↓'}},
    '지질억제': {'genes': ['SREBPA', 'SREBPa', 'SREBPC', 'SREBPc'], 'cell': 'SZ95',
              'expected': {'SREBPA': '↓', 'SREBPC': '↓'}},
    '냉감': {'genes': ['TRPM8', 'CIRBP'], 'cell': 'HaCaT'}
}

# Valid colors
VALID_COLORS = {
    'black': '#000000', 'white': '#FFFFFF', 'gray': '#808080', 'lightgray': '#D3D3D3',
    'red': '#FF0000', 'blue': '#0000FF', 'green': '#008000', 'yellow': '#FFFF00'
}

# ==================== UNIVERSAL PARSER ====================
class UniversalParser:
    @staticmethod
    def similarity(a, b):
        """Calculate string similarity"""
        return SequenceMatcher(None, str(a).lower(), str(b).lower()).ratio()
    
    @staticmethod
    def find_data_start(df):
        """Find where actual data starts by looking for 'Well' or 'Sample' keywords"""
        keywords = ['well', 'sample', 'target', 'ct', 'cт']
        
        for idx in range(min(50, len(df))):
            row_str = ' '.join(str(x).lower() for x in df.iloc[idx].values if pd.notna(x))
            if any(kw in row_str for kw in keywords):
                return idx
        return 0
    
    @staticmethod
    def find_columns(df):
        """Intelligently find Well, Sample, Target, CT columns"""
        header_row = df.iloc[0]
        
        well_col = None
        sample_col = None
        target_col = None
        ct_col = None
        
        # Find columns by fuzzy matching
        for i, col_name in enumerate(header_row):
            col_lower = str(col_name).lower()
            
            if well_col is None and ('well' in col_lower):
                well_col = i
            elif sample_col is None and ('sample' in col_lower):
                sample_col = i
            elif target_col is None and ('target' in col_lower or 'gene' in col_lower):
                target_col = i
            elif ct_col is None and ('ct' in col_lower or 'cт' in col_lower or 'cq' in col_lower):
                ct_col = i
        
        return well_col, sample_col, target_col, ct_col
    
    @staticmethod
    def parse(file):
        """Universal CSV parser - works with any format"""
        try:
            # Try multiple encodings
            df = None
            for enc in ['utf-8', 'latin-1', 'cp1252', 'iso-8859-1']:
                try:
                    df = pd.read_csv(file, encoding=enc, header=None, skip_blank_lines=False)
                    break
                except:
                    continue
            
            if df is None or len(df) == 0:
                return None
            
            # Find data start
            start_idx = UniversalParser.find_data_start(df)
            df = df.iloc[start_idx:].reset_index(drop=True)
            
            # Set first row as header
            df.columns = df.iloc[0]
            df = df.iloc[1:].reset_index(drop=True)
            
            # Find columns
            well_idx, sample_idx, target_idx, ct_idx = UniversalParser.find_columns(df)
            
            # Fallback to position if not found
            if sample_idx is None:
                sample_idx = 2 if len(df.columns) > 2 else 1
            if target_idx is None:
                target_idx = 3 if len(df.columns) > 3 else 2
            if ct_idx is None:
                # Find first numeric column
                for i in range(len(df.columns)):
                    if pd.to_numeric(df.iloc[:, i], errors='coerce').notna().sum() > len(df) * 0.5:
                        ct_idx = i
                        break
            
            if ct_idx is None:
                return None
            
            # Extract data
            parsed = pd.DataFrame({
                'Well': df.iloc[:, well_idx] if well_idx is not None else df.iloc[:, 0],
                'Sample': df.iloc[:, sample_idx],
                'Target': df.iloc[:, target_idx],
                'CT': pd.to_numeric(df.iloc[:, ct_idx], errors='coerce')
            })
            
            # Clean data
            parsed = parsed.dropna(subset=['CT'])
            parsed = parsed[parsed['Sample'].notna() & parsed['Target'].notna()]
            parsed = parsed[parsed['Sample'].astype(str).str.strip() != '']
            
            return parsed if len(parsed) > 0 else None
            
        except Exception as e:
            st.error(f"Parse error: {str(e)}")
            return None
    
    @staticmethod
    def auto_merge_samples(df, threshold=0.8):
        """Automatically detect and suggest sample merges"""
        samples = df['Sample'].unique()
        merge_suggestions = {}
        
        for i, s1 in enumerate(samples):
            for s2 in samples[i+1:]:
                sim = UniversalParser.similarity(s1, s2)
                if sim > threshold:
                    # Use shorter name or first one
                    canonical = s1 if len(str(s1)) <= len(str(s2)) else s2
                    merge_suggestions[s1] = canonical
                    merge_suggestions[s2] = canonical
        
        return merge_suggestions

# ==================== ANALYSIS ENGINE ====================
class AnalysisEngine:
    @staticmethod
    def calculate_ddct(data, hk_gene, baseline, mapping, included):
        data = data.copy()
        data['Condition'] = data['Sample'].map(lambda x: mapping.get(x, {}).get('condition', x))
        data['Group'] = data['Sample'].map(lambda x: mapping.get(x, {}).get('group', 'Treatment'))
        
        results = {}
        
        for target in data['Target'].unique():
            if target.upper() in [hk_gene.upper(), 'ACTIN', 'B-ACTIN', 'GAPDH', 'ACTB']:
                continue
            
            inc = included.get(target, [])
            if not inc or baseline not in inc:
                continue
            
            tdata = data[(data['Target'] == target) & (data['Condition'].isin(inc))]
            
            # Baseline
            bl_t = tdata[tdata['Condition'] == baseline]
            bl_h = data[(data['Condition'] == baseline) & (data['Target'] == hk_gene)]
            if len(bl_t) == 0 or len(bl_h) == 0:
                continue
            
            bl_dct = bl_t['CT'].mean() - bl_h['CT'].mean()
            
            gene_results = []
            for cond in tdata['Condition'].unique():
                ct = tdata[tdata['Condition'] == cond]
                ch = data[(data['Condition'] == cond) & (data['Target'] == hk_gene)]
                
                if len(ch) == 0:
                    continue
                
                ct_vals = ct['CT'].values
                dct = ct_vals.mean() - ch['CT'].mean()
                ddct = dct - bl_dct
                fc = 2 ** (-ddct)
                sem = ct_vals.std() / np.sqrt(len(ct_vals)) if len(ct_vals) > 1 else 0
                
                gene_results.append({
                    'Target': target, 'Condition': cond, 'Group': ct['Group'].iloc[0],
                    'n': len(ct_vals), 'CT_values': ct_vals, 'Fold_Change': fc, 'SEM': sem
                })
            
            if gene_results:
                results[target] = pd.DataFrame(gene_results)
        
        return results
    
    @staticmethod
    def calculate_stats(df, ref):
        df = df.copy()
        df['p_value'] = np.nan
        df['sig'] = ''
        
        ref_data = df[df['Condition'] == ref]
        if len(ref_data) == 0:
            return df
        
        ref_vals = ref_data.iloc[0]['CT_values']
        
        for idx, row in df.iterrows():
            if row['Condition'] == ref:
                df.at[idx, 'p_value'] = 1.0
                continue
            
            vals = row['CT_values']
            
            try:
                if len(ref_vals) > 1 and len(vals) > 1:
                    _, p = stats.ttest_ind(ref_vals, vals)
                elif len(ref_vals) == 1 and len(vals) > 1:
                    _, p = stats.ttest_1samp(vals, ref_vals[0])
                else:
                    p = np.nan
                
                df.at[idx, 'p_value'] = p
                
                if p < 0.001:
                    df.at[idx, 'sig'] = '***'
                elif p < 0.01:
                    df.at[idx, 'sig'] = '**'
                elif p < 0.05:
                    df.at[idx, 'sig'] = '*'
            except:
                pass
        
        return df

# ==================== GRAPH GENERATOR ====================
def create_graph(data, gene, settings):
    gdata = data[data['Condition'].isin(settings['show_conditions'])].copy()
    
    # Custom order
    if settings.get('order'):
        order_map = {c: i for i, c in enumerate(settings['order'])}
        gdata['_sort'] = gdata['Condition'].map(lambda x: order_map.get(x, 999))
        gdata = gdata.sort_values('_sort')
    
    colors = [settings['colors'].get(c, '#636EFA') for c in gdata['Condition']]
    
    fig = go.Figure()
    
    fig.add_trace(go.Bar(
        x=gdata['Condition'],
        y=gdata['Fold_Change'],
        error_y=dict(
            type='data',
            array=gdata['SEM'] * 1.96,
            visible=settings.get('show_error', True),
            thickness=settings.get('err_thick', 2),
            width=settings.get('err_width', 4),
            color=settings.get('err_color', '#000000')
        ),
        text=gdata['sig'] if settings.get('show_sig', True) else None,
        textposition=settings.get('sig_pos', 'outside'),
        textfont=dict(size=settings.get('sig_size', 16), color=settings.get('sig_color', '#000000')),
        marker=dict(
            color=colors,
            line=dict(color=settings.get('border_color', '#000000'), width=settings.get('border_width', 0)),
            opacity=settings.get('opacity', 1.0)
        ),
        width=settings.get('bar_width', 0.8),
        showlegend=False
    ))
    
    # Reference line
    if settings.get('show_ref', True):
        fig.add_hline(y=settings.get('ref_val', 1.0), line_dash="dash",
                     line_color=settings.get('ref_color', '#808080'),
                     line_width=settings.get('ref_width', 2))
    
    fig.update_layout(
        title=dict(text=settings.get('title', gene), font=dict(size=settings.get('title_size', 20))),
        xaxis=dict(
            title=settings.get('xlabel', 'Condition'),
            showgrid=settings.get('x_grid', False),
            gridcolor=VALID_COLORS.get('lightgray', '#D3D3D3'),
            tickangle=settings.get('tick_angle', -45),
            tickfont=dict(size=settings.get('tick_size', 12))
        ),
        yaxis=dict(
            title=settings.get('ylabel', 'Fold Change'),
            showgrid=settings.get('y_grid', True),
            gridcolor=VALID_COLORS.get('lightgray', '#D3D3D3'),
            range=settings.get('y_range')
        ),
        template=settings.get('theme', 'plotly_white'),
        height=settings.get('height', 600),
        width=settings.get('width', 1000),
        plot_bgcolor=settings.get('plot_bg', '#FFFFFF'),
        paper_bgcolor=settings.get('paper_bg', '#FFFFFF')
    )
    
    return fig

# ==================== EXPORT ====================
def export_excel(raw, processed, params, mapping):
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine='xlsxwriter') as w:
        pd.DataFrame([params]).to_excel(w, sheet_name='Params', index=False)
        pd.DataFrame([{'Original': k, **v} for k, v in mapping.items()]).to_excel(w, sheet_name='Mapping', index=False)
        raw.to_excel(w, sheet_name='Raw', index=False)
        for g, d in processed.items():
            d.drop(columns=['CT_values'], errors='ignore').to_excel(w, sheet_name=g[:31], index=False)
    return output.getvalue()

# ==================== UI ====================
st.title("🧬 qPCR Analysis Ultimate Pro")

with st.sidebar:
    st.header("📖 Quick Guide")
    st.markdown("""
    1. Upload CSV files
    2. Auto-merge samples
    3. Map conditions
    4. Analyze genes
    5. Customize graphs
    6. Export results
    """)
    
    if st.button("💾 Save Template"):
        name = st.text_input("Name:", key="tpl_name")
        if name and st.session_state.sample_mapping:
            st.session_state.analysis_templates[name] = {
                'mapping': st.session_state.sample_mapping.copy(),
                'settings': st.session_state.graph_settings.copy()
            }
            st.success(f"✅ Saved '{name}'")

tabs = st.tabs(["📁 Data", "🗺️ Mapping", "🔬 Analysis", "📊 Graphs", "📤 Export"])

# ==================== TAB 1: DATA ====================
with tabs[0]:
    st.header("Upload & Process Data")
    
    files = st.file_uploader("Upload CSV files", type=['csv'], accept_multiple_files=True)
    
    if files:
        all_data = []
        for f in files:
            parsed = UniversalParser.parse(f)
            if parsed is not None:
                parsed['Source'] = f.name
                all_data.append(parsed)
                st.success(f"✅ {f.name}: {len(parsed)} wells")
        
        if all_data:
            combined = pd.concat(all_data, ignore_index=True)
            
            # Auto-merge detection
            st.subheader("🔗 Auto-Merge Similar Samples")
            merge_sugg = UniversalParser.auto_merge_samples(combined)
            
            if merge_sugg:
                st.info(f"Found {len(set(merge_sugg.values()))} potential merges")
                
                if st.checkbox("Review merge suggestions"):
                    for orig, canonical in sorted(set(merge_sugg.items())):
                        col1, col2, col3 = st.columns([2, 1, 2])
                        col1.text(f"'{orig}'")
                        col2.markdown("→")
                        apply = col3.checkbox(f"Merge to '{canonical}'", key=f"m_{orig}")
                        if apply:
                            st.session_state.auto_merged[orig] = canonical
                
                if st.button("Apply Auto-Merges") and st.session_state.auto_merged:
                    combined['Sample'] = combined['Sample'].map(
                        lambda x: st.session_state.auto_merged.get(x, x)
                    )
                    st.success(f"✅ Merged {len(st.session_state.auto_merged)} samples")
            
            st.session_state.data = combined
            
            # HK gene
            hk = [g for g in combined['Target'].unique() if g.upper() in ['ACTIN', 'B-ACTIN', 'GAPDH', 'ACTB']]
            if hk:
                st.session_state.hk_gene = st.selectbox("🔬 Housekeeping Gene", hk)
            
            # Filters
            st.subheader("🎯 Filters")
            samples = sorted(combined['Sample'].unique())
            sel_samples = st.multiselect("Samples:", samples, default=samples, key="filt_samples")
            st.session_state.excluded_samples = set(samples) - set(sel_samples)
            
            genes = [g for g in sorted(combined['Target'].unique()) if g.upper() not in ['ACTIN', 'GAPDH', 'B-ACTIN', 'ACTB']]
            sel_genes = st.multiselect("Genes:", genes, default=genes, key="filt_genes")
            
            # Show data
            st.subheader("📊 Filtered Data")
            filt = combined[combined['Sample'].isin(sel_samples) & combined['Target'].isin(sel_genes)]
            
            for s in sel_samples[:5]:  # Show first 5
                with st.expander(f"📍 {s}"):
                    sdata = filt[filt['Sample'] == s]
                    st.dataframe(sdata.groupby('Target')['CT'].agg(['mean', 'std', 'count']).round(2))

# ==================== TAB 2: MAPPING ====================
with tabs[1]:
    st.header("Sample Mapping")
    
    if st.session_state.data is not None:
        detected = set(st.session_state.data['Target'].unique())
        suggested = next((e for e, c in EFFICACY_CONFIG.items() if any(g in detected for g in c['genes'])), None)
        
        efficacy = st.selectbox("🎯 Efficacy Type", list(EFFICACY_CONFIG.keys()),
                               index=list(EFFICACY_CONFIG.keys()).index(suggested) if suggested else 0)
        st.session_state.selected_efficacy = efficacy
        
        st.subheader("🗺️ Map Conditions")
        
        samples = [s for s in sorted(st.session_state.data['Sample'].unique()) 
                  if s not in st.session_state.excluded_samples]
        
        for s in samples:
            if s not in st.session_state.sample_mapping:
                st.session_state.sample_mapping[s] = {'condition': s, 'group': 'Treatment', 'conc': ''}
        
        # Headers
        c1, c2, c3, c4 = st.columns([1, 3, 2, 2])
        c1.markdown("**Original**")
        c2.markdown("**Condition**")
        c3.markdown("**Group**")
        c4.markdown("**Concentration**")
        
        groups = ['Baseline', 'Negative Control', 'Positive Control', 'Treatment']
        
        for s in samples:
            c1, c2, c3, c4 = st.columns([1, 3, 2, 2])
            c1.text(s)
            
            cond = c2.text_input("C", st.session_state.sample_mapping[s]['condition'],
                                key=f"cn_{s}", label_visibility="collapsed")
            st.session_state.sample_mapping[s]['condition'] = cond
            
            grp = c3.selectbox("G", groups,
                              index=groups.index(st.session_state.sample_mapping[s]['group']) 
                              if st.session_state.sample_mapping[s]['group'] in groups else 0,
                              key=f"gr_{s}", label_visibility="collapsed")
            st.session_state.sample_mapping[s]['group'] = grp
            
            conc = c4.text_input("Cn", st.session_state.sample_mapping[s]['conc'],
                                key=f"cc_{s}", label_visibility="collapsed")
            st.session_state.sample_mapping[s]['conc'] = conc

# ==================== TAB 3: ANALYSIS ====================
with tabs[2]:
    st.header("Analysis")
    
    if st.session_state.data is not None and st.session_state.sample_mapping and st.session_state.hk_gene:
        c1, c2 = st.columns(2)
        
        with c1:
            bl_samples = [k for k, v in st.session_state.sample_mapping.items() if v['group'] == 'Baseline']
            bl = st.selectbox("🎯 Baseline", bl_samples if bl_samples else list(st.session_state.sample_mapping.keys()))
            bl_cond = st.session_state.sample_mapping[bl]['condition']
        
        with c2:
            all_conds = [v['condition'] for v in st.session_state.sample_mapping.values()]
            pv_cond = st.selectbox("📊 P-value Ref", all_conds,
                                  index=all_conds.index(bl_cond) if bl_cond in all_conds else 0)
        
        # Per-gene selection
        st.subheader("🧬 Samples Per Gene")
        genes = [g for g in st.session_state.data['Target'].unique() 
                if g.upper() not in ['ACTIN', 'GAPDH', 'B-ACTIN', 'ACTB']]
        
        for g in genes:
            if g not in st.session_state.included_samples_analysis:
                st.session_state.included_samples_analysis[g] = all_conds.copy()
            
            sel = st.multiselect(f"{g}:", all_conds,
                                default=st.session_state.included_samples_analysis[g],
                                key=f"inc_{g}")
            st.session_state.included_samples_analysis[g] = sel
        
        if st.button("🔬 Run Analysis", type="primary"):
            with st.spinner("Analyzing..."):
                results = AnalysisEngine.calculate_ddct(
                    st.session_state.data, st.session_state.hk_gene, bl_cond,
                    st.session_state.sample_mapping, st.session_state.included_samples_analysis
                )
                
                if results:
                    for g, d in results.items():
                        results[g] = AnalysisEngine.calculate_stats(d, pv_cond)
                    
                    st.session_state.processed_data = results
                    st.success(f"✅ {len(results)} genes analyzed")
        
        if st.session_state.processed_data:
            st.subheader("📊 Results")
            for g, d in st.session_state.processed_data.items():
                with st.expander(f"{g}", expanded=False):
                    disp = d.drop(columns=['CT_values'], errors='ignore')
                    st.dataframe(disp[['Condition', 'Group', 'Fold_Change', 'p_value', 'sig', 'n']]\
                        .style.background_gradient(subset=['Fold_Change'], cmap='RdYlGn')\
                        .format({'Fold_Change': '{:.3f}', 'p_value': '{:.4f}'}))

# ==================== TAB 4: GRAPHS ====================
with tabs[3]:
    st.header("Interactive Graph Editor")
    
    if st.session_state.processed_data:
        
        # Quick presets
        c1, c2, c3 = st.columns(3)
        if c1.button("📄 Publication"):
            for g in st.session_state.graph_settings:
                st.session_state.graph_settings[g].update({'theme': 'simple_white', 'width': 1000, 'height': 600})
            st.rerun()
        if c2.button("🎨 Presentation"):
            for g in st.session_state.graph_settings:
                st.session_state.graph_settings[g].update({'theme': 'presentation', 'width': 1200, 'height': 700})
            st.rerun()
        if c3.button("🔄 Reset All"):
            st.session_state.graph_settings = {}
            st.rerun()
        
        st.markdown("---")
        
        for gene in st.session_state.processed_data.keys():
            st.markdown(f"## {gene}")
            
            # Initialize
            if gene not in st.session_state.graph_settings:
                colors_default = px.colors.qualitative.Plotly
                all_conds = list(st.session_state.processed_data[gene]['Condition'].unique())
                
                st.session_state.graph_settings[gene] = {
                    'title': f"{gene} Expression", 'xlabel': 'Condition', 'ylabel': 'Fold Change',
                    'title_size': 20, 'tick_size': 12, 'sig_size': 16,
                    'width': 1000, 'height': 600, 'theme': 'plotly_white',
                    'show_error': True, 'show_sig': True, 'show_ref': True, 'ref_val': 1.0,
                    'x_grid': False, 'y_grid': True, 'bar_width': 0.8, 'opacity': 1.0,
                    'border_width': 0, 'err_thick': 2, 'err_width': 4, 'ref_width': 2,
                    'tick_angle': -45, 'sig_pos': 'outside', 'y_range': None,
                    'show_conditions': all_conds, 'order': None,
                    'colors': {c: colors_default[i % len(colors_default)] for i, c in enumerate(all_conds)},
                    'sig_color': '#000000', 'err_color': '#000000', 'border_color': '#000000',
                    'ref_color': '#808080', 'plot_bg': '#FFFFFF', 'paper_bg': '#FFFFFF'
                }
            
            settings = st.session_state.graph_settings[gene]
            gdata = st.session_state.processed_data[gene]
            
            # Tabs
            t1, t2, t3, t4 = st.tabs(["📋 Data", "📊 Style", "🎨 Colors", "📝 Text"])
            
            with t1:
                all_c = list(gdata['Condition'].unique())
                settings['show_conditions'] = st.multiselect("Display:", all_c, default=settings['show_conditions'],
                                                             key=f"show_{gene}")
                
                if st.checkbox("Custom order", key=f"ord_{gene}"):
                    settings['order'] = st.multiselect("Order:", settings['show_conditions'],
                                                      default=settings['show_conditions'], key=f"order_{gene}")
            
            with t2:
                c1, c2, c3 = st.columns(3)
                with c1:
                    settings['bar_width'] = st.slider("Bar Width", 0.1, 1.5, float(settings['bar_width']), key=f"bw_{gene}")
                    settings['opacity'] = st.slider("Opacity", 0.0, 1.0, float(settings['opacity']), key=f"op_{gene}")
                    settings['border_width'] = st.slider("Border", 0, 5, int(settings['border_width']), key=f"bd_{gene}")
                
                with c2:
                    settings['show_error'] = st.checkbox("Error Bars", settings['show_error'], key=f"er_{gene}")
                    if settings['show_error']:
                        settings['err_thick'] = st.slider("Thickness", 1, 5, int(settings['err_thick']), key=f"et_{gene}")
                        settings['err_width'] = st.slider("Cap Width", 2, 10, int(settings['err_width']), key=f"ew_{gene}")
                
                with c3:
                    settings['show_ref'] = st.checkbox("Reference Line", settings['show_ref'], key=f"rf_{gene}")
                    if settings['show_ref']:
                        settings['ref_val'] = st.number_input("Value", value=float(settings['ref_val']), key=f"rv_{gene}")
                        settings['ref_width'] = st.slider("Width", 1, 5, int(settings['ref_width']), key=f"rw_{gene}")
                
                c1, c2, c3 = st.columns(3)
                with c1:
                    themes = ['plotly_white', 'plotly', 'plotly_dark', 'seaborn', 'simple_white', 'presentation']
                    settings['theme'] = st.selectbox("Theme", themes, 
                                                     index=themes.index(settings['theme']) if settings['theme'] in themes else 0,
                                                     key=f"th_{gene}")
                with c2:
                    settings['width'] = st.number_input("Width", 400, 2000, int(settings['width']), key=f"w_{gene}")
                    settings['height'] = st.number_input("Height", 300, 1500, int(settings['height']), key=f"h_{gene}")
                with c3:
                    settings['x_grid'] = st.checkbox("X Grid", settings['x_grid'], key=f"xg_{gene}")
                    settings['y_grid'] = st.checkbox("Y Grid", settings['y_grid'], key=f"yg_{gene}")
                    
                    use_y = st.checkbox("Custom Y Range", value=settings['y_range'] is not None, key=f"yr_{gene}")
                    if use_y:
                        curr = settings['y_range'] or [0, 3]
                        ymin = st.number_input("Y Min", value=float(curr[0]), key=f"ymin_{gene}")
                        ymax = st.number_input("Y Max", value=float(curr[1]), key=f"ymax_{gene}")
                        settings['y_range'] = [ymin, ymax]
                    else:
                        settings['y_range'] = None
            
            with t3:
                st.markdown("**Bar Colors**")
                cols = st.columns(3)
                for i, c in enumerate(settings['show_conditions']):
                    with cols[i % 3]:
                        if c not in settings['colors']:
                            settings['colors'][c] = '#636EFA'
                        settings['colors'][c] = st.color_picker(c, settings['colors'][c], key=f"col_{gene}_{c}")
                
                st.markdown("**Other Colors**")
                c1, c2, c3 = st.columns(3)
                with c1:
                    settings['sig_color'] = st.color_picker("Significance", settings['sig_color'], key=f"sc_{gene}")
                    settings['err_color'] = st.color_picker("Error Bars", settings['err_color'], key=f"ec_{gene}")
                with c2:
                    settings['border_color'] = st.color_picker("Border", settings['border_color'], key=f"bc_{gene}")
                    settings['ref_color'] = st.color_picker("Reference", settings['ref_color'], key=f"rc_{gene}")
                with c3:
                    settings['plot_bg'] = st.color_picker("Plot BG", settings['plot_bg'], key=f"pb_{gene}")
                    settings['paper_bg'] = st.color_picker("Paper BG", settings['paper_bg'], key=f"papb_{gene}")
            
            with t4:
                c1, c2 = st.columns(2)
                with c1:
                    settings['title'] = st.text_input("Title", settings['title'], key=f"tit_{gene}")
                    settings['xlabel'] = st.text_input("X Label", settings['xlabel'], key=f"xl_{gene}")
                    settings['ylabel'] = st.text_input("Y Label", settings['ylabel'], key=f"yl_{gene}")
                
                with c2:
                    settings['title_size'] = st.slider("Title Size", 12, 36, int(settings['title_size']), key=f"ts_{gene}")
                    settings['tick_size'] = st.slider("Tick Size", 8, 20, int(settings['tick_size']), key=f"tks_{gene}")
                    settings['tick_angle'] = st.slider("Tick Angle", -90, 90, int(settings['tick_angle']), key=f"ta_{gene}")
                
                c1, c2 = st.columns(2)
                with c1:
                    settings['show_sig'] = st.checkbox("Show Stars", settings['show_sig'], key=f"ss_{gene}")
                with c2:
                    if settings['show_sig']:
                        positions = ['outside', 'inside', 'auto']
                        settings['sig_pos'] = st.selectbox("Position", positions,
                                                          index=positions.index(settings['sig_pos']) if settings['sig_pos'] in positions else 0,
                                                          key=f"sp_{gene}")
                        settings['sig_size'] = st.slider("Star Size", 10, 30, int(settings['sig_size']), key=f"ssi_{gene}")
            
            # Generate graph
            st.markdown("### 📊 Preview")
            fig = create_graph(gdata, gene, settings)
            st.plotly_chart(fig, use_container_width=True, key=f"fig_{gene}")
            
            if 'graphs' not in st.session_state:
                st.session_state.graphs = {}
            st.session_state.graphs[gene] = fig
            
            st.markdown("---")

# ==================== TAB 5: EXPORT ====================
with tabs[4]:
    st.header("Export Results")
    
    if st.session_state.processed_data:
        params = {
            'Date': datetime.now().strftime("%Y-%m-%d %H:%M"),
            'Efficacy': st.session_state.selected_efficacy,
            'HK_Gene': st.session_state.hk_gene,
            'Baseline': bl_cond if 'bl_cond' in locals() else 'N/A',
            'Pvalue_Ref': pv_cond if 'pv_cond' in locals() else 'N/A',
            'Genes': len(st.session_state.processed_data)
        }
        
        c1, c2 = st.columns(2)
        
        with c1:
            st.markdown("### 📊 Excel")
            excel = export_excel(st.session_state.data, st.session_state.processed_data,
                                params, st.session_state.sample_mapping)
            st.download_button("📥 Download Excel", excel,
                             f"qPCR_{datetime.now().strftime('%Y%m%d')}.xlsx",
                             "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                             type="primary")
        
        with c2:
            st.markdown("### 📈 Graphs")
            if st.session_state.graphs:
                html = ["<html><head><title>qPCR</title></head><body>",
                       f"<h1>{st.session_state.selected_efficacy}</h1>"]
                for g, fig in st.session_state.graphs.items():
                    html.append(f"<h2>{g}</h2>")
                    html.append(fig.to_html(include_plotlyjs='cdn'))
                html.append("</body></html>")
                
                st.download_button("📥 Download Graphs", "\n".join(html),
                                 f"graphs_{datetime.now().strftime('%Y%m%d')}.html",
                                 "text/html", type="primary")
        
        st.markdown("---")
        st.subheader("Individual Files")
        
        c1, c2, c3 = st.columns(3)
        
        with c1:
            st.markdown("**CSV**")
            for g, d in st.session_state.processed_data.items():
                csv = io.StringIO()
                d.drop(columns=['CT_values'], errors='ignore').to_csv(csv, index=False)
                st.download_button(f"{g}.csv", csv.getvalue(), f"{g}.csv", "text/csv", key=f"csv_{g}")
        
        with c2:
            st.markdown("**HTML**")
            for g, fig in st.session_state.graphs.items():
                html = io.StringIO()
                fig.write_html(html)
                st.download_button(f"{g}.html", html.getvalue(), f"{g}.html", "text/html", key=f"html_{g}")
        
        with c3:
            st.markdown("**Config**")
            config = {'params': params, 'mapping': st.session_state.sample_mapping,
                     'settings': st.session_state.graph_settings}
            st.download_button("Config", json.dumps(config, indent=2), "config.json", "application/json")

st.markdown("---")
st.caption("🧬 qPCR Analysis Ultimate Pro v5.0 | Universal Parser • Auto-Merge • Full Customization")
