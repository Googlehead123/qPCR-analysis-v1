import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from scipy import stats
import io
import warnings
import json
import re
import zipfile
from datetime import datetime
from typing import Dict, Tuple


# ==================== UTILITY FUNCTIONS ====================
def natural_sort_key(sample_name):
    """Extract numbers from sample name for natural sorting (e.g., Sample2 < Sample10)"""
    parts = re.split(r"(\d+)", str(sample_name))
    return [int(part) if part.isdigit() else part.lower() for part in parts]


# ==================== PAGE CONFIG ====================
st.set_page_config(
    page_title="qPCR Analysis Suite Pro",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ==================== CONSTANTS ====================
DEFAULT_GROUP_COLORS = {
    "Baseline": "#FFFFFF",
    "Non-treated": "#FFFFFF",
    "Control": "#FFFFFF",
    "Negative Control": "#9E9E9E",
    "Inducer": "#9E9E9E",
    "Positive Control": "#9E9E9E",
    "Treatment": "#D3D3D3",
}

# ==================== SESSION STATE INIT ====================
for key in [
    "data",
    "processed_data",
    "sample_mapping",
    "analysis_templates",
    "graphs",
    "excluded_wells",
    "excluded_samples",
    "selected_efficacy",
    "hk_gene",
]:
    if key not in st.session_state:
        st.session_state[key] = (
            {}
            if key in ["sample_mapping", "analysis_templates", "graphs"]
            else (set() if "excluded" in key else None)
        )

# ==================== EFFICACY DATABASE ====================
EFFICACY_CONFIG = {
    "탄력": {
        "genes": ["COL1A1", "ELN", "FBN-1", "FBN1"],
        "cell": "HS68 fibroblast",
        "controls": {
            "negative": "Non-treated",
            "positive": "TGFb",
            "compare_to": "negative",
        },
        "description": "Elasticity - Non-treated vs TGFb (positive) vs Treatments",
    },
    "항노화": {
        "genes": ["COL1A1", "COL1", "MMP-1", "MMP1"],
        "cell": "HS68 fibroblast",
        "controls": {
            "baseline": "Non-treated (No UV)",
            "negative": "UVB only",
            "positive": "UVB+TGFb",
            "compare_to": "negative",
        },
        "description": "Anti-aging - COL1↑ (recovery), MMP1↓ (inhibition) after UVB damage",
        "expected_direction": {
            "COL1A1": "up",
            "COL1": "up",
            "MMP-1": "down",
            "MMP1": "down",
        },
    },
    "보습": {
        "genes": ["AQP3", "HAS3"],
        "cell": "HaCaT keratinocyte",
        "controls": {
            "negative": "Non-treated",
            "positive": "Retinoic acid",
            "compare_to": "negative",
        },
        "description": "Hydration - Non-treated vs Retinoic acid (positive) vs Treatments",
    },
    "장벽": {
        "genes": ["FLG", "CLDN", "IVL"],
        "cell": "HaCaT keratinocyte",
        "controls": {
            "negative": "Non-treated",
            "positive": "Retinoic acid",
            "compare_to": "negative",
        },
        "description": "Barrier function - Non-treated vs Retinoic acid (positive) vs Treatments",
    },
    "표피증식": {
        "genes": ["KI67", "PCNA"],
        "cell": "HaCaT keratinocyte",
        "controls": {
            "negative": "Non-treated",
            "positive": "TGFb or FBS",
            "compare_to": "negative",
        },
        "description": "Proliferation - Non-treated vs TGFb/FBS (positive) vs Treatments",
    },
    "멜라닌억제": {
        "genes": ["MITF", "TYR", "Melanin"],
        "cell": "B16F10 melanocyte",
        "controls": {
            "baseline": "Non-treated",
            "negative": "α-MSH only",
            "positive": "α-MSH+Arbutin",
            "compare_to": "negative",
        },
        "description": "Melanin inhibition - α-MSH induced vs α-MSH+Arbutin (positive) vs α-MSH+Treatments",
        "expected_direction": {"MITF": "down", "TYR": "down", "Melanin": "down"},
    },
    "진정": {
        "genes": ["IL1B", "IL-1β", "IL6", "TNFA", "TNFα"],
        "cell": "HaCaT keratinocyte",
        "controls": {
            "baseline": "Non-treated",
            "negative": "IL4+PolyIC (Inflammation)",
            "positive": "Inflammation+Dexamethasone",
            "compare_to": "negative",
        },
        "description": "Anti-inflammation - Reduce IL1β/IL6/TNFα (all should decrease)",
        "expected_direction": {
            "IL1B": "down",
            "IL-1β": "down",
            "IL6": "down",
            "TNFA": "down",
            "TNFα": "down",
        },
    },
    "지질억제": {
        "genes": ["SREBPA", "SREBPa", "SREBPC", "SREBPc", "PPARY", "PPARy"],
        "cell": "SZ95 sebocyte",
        "controls": {
            "baseline": "Non-treated",
            "negative": "IGF only",
            "positive": "IGF+Reference inhibitor",
            "compare_to": "negative",
        },
        "description": "Sebum inhibition - IGF induced vs IGF+Treatments",
        "expected_direction": {
            "SREBPA": "down",
            "SREBPa": "down",
            "SREBPC": "down",
            "SREBPc": "down",
            "PPARY": "down",
            "PPARy": "down",
        },
    },
    "냉감": {
        "genes": ["TRPM8", "CIRBP"],
        "cell": "HaCaT keratinocyte",
        "controls": {
            "negative": "Non-treated",
            "positive": "Menthol",
            "compare_to": "negative",
        },
        "description": "Cooling effect - Non-treated vs Menthol (positive) vs Treatments",
    },
}


# ==================== ANALYSIS CONSTANTS ====================
class AnalysisConstants:
    MIN_REPLICATES_FOR_STATS = 2
    RECOMMENDED_REPLICATES = 3
    CT_UNDETERMINED_THRESHOLD = 40.0
    CT_HIGH_WARNING = 35.0
    CT_LOW_WARNING = 10.0
    CV_WARNING_THRESHOLD = 0.05
    CV_ERROR_THRESHOLD = 0.10
    HK_DEVIATION_WARNING = 1.0
    HK_DEVIATION_ERROR = 2.0
    P_VALUE_THRESHOLDS = {"*": 0.05, "**": 0.01, "***": 0.001}


# ==================== PARSER CLASS ====================
class QPCRParser:
    @staticmethod
    def detect_format(df):
        for idx, row in df.iterrows():
            row_str = " ".join(row.astype(str).values)
            if "Well Position" in row_str:
                return "format1", idx
            elif row.iloc[0] == "Well" and "Sample Name" in row_str:
                return (
                    "format2" if "Cт" in row_str or "ΔCт" in row_str else "format1"
                ), idx
        return "unknown", 0

    @staticmethod
    def parse_format1(df, start):
        df = df.iloc[start:].reset_index(drop=True)
        df.columns = df.iloc[0]
        df = df.iloc[1:].reset_index(drop=True)

        if len(df.columns) < 4:
            st.error("CSV must have at least 4 columns (Well, Sample, Target, CT)")
            return None

        well_col = next(
            (c for c in ["Well Position", "Well"] if c in df.columns), df.columns[0]
        )
        ct_col = next((c for c in ["CT", "Ct", "Cт"] if c in df.columns), None)

        if not ct_col:
            return None

        sample_col = (
            df.get("Sample Name", df.iloc[:, 2])
            if "Sample Name" in df.columns
            else df.iloc[:, 2]
        )
        target_col = (
            df.get("Target Name", df.iloc[:, 3])
            if "Target Name" in df.columns
            else df.iloc[:, 3]
        )

        parsed = pd.DataFrame(
            {
                "Well": df[well_col],
                "Sample": sample_col,
                "Target": target_col,
                "CT": pd.to_numeric(df[ct_col], errors="coerce"),
            }
        )

        return parsed.dropna(subset=["CT"]).query("Sample.notna() & Target.notna()")

    @staticmethod
    def parse_format2(df, start):
        df = df.iloc[start:].reset_index(drop=True)
        df.columns = df.iloc[0]
        df = df.iloc[1:].reset_index(drop=True)

        parsed = pd.DataFrame(
            {
                "Well": df["Well"],
                "Sample": df["Sample Name"],
                "Target": df["Target Name"],
                "CT": pd.to_numeric(df["Cт"], errors="coerce"),
            }
        )

        return parsed.dropna(subset=["CT"]).query("Sample.notna() & Target.notna()")

    MAX_FILE_SIZE_MB = 50

    @staticmethod
    def parse(file):
        try:
            file.seek(0, 2)
            file_size_mb = file.tell() / (1024 * 1024)
            file.seek(0)

            if file_size_mb > QPCRParser.MAX_FILE_SIZE_MB:
                st.error(
                    f"File too large ({file_size_mb:.1f} MB). Maximum size is {QPCRParser.MAX_FILE_SIZE_MB} MB."
                )
                return None

            df = None
            for enc in ["utf-8", "latin-1", "cp1252"]:
                try:
                    df = pd.read_csv(
                        file, encoding=enc, low_memory=False, skip_blank_lines=False
                    )
                    break
                except UnicodeDecodeError:
                    file.seek(0)
                    continue

            if df is None:
                return None

            fmt, start = QPCRParser.detect_format(df)
            return (
                QPCRParser.parse_format1(df, start)
                if fmt == "format1"
                else QPCRParser.parse_format2(df, start)
                if fmt == "format2"
                else None
            )
        except Exception as e:
            st.error(f"Parse error: {e}")
            return None


# ==================== QUALITY CONTROL ====================
class QualityControl:
    CT_HIGH_THRESHOLD = 35.0
    CT_LOW_THRESHOLD = 10.0
    CV_THRESHOLD = 0.05
    HK_VARIATION_THRESHOLD = 1.0
    GRUBBS_ALPHA = 0.05

    @staticmethod
    def get_triplicate_data(
        data: pd.DataFrame, excluded_wells: set = None
    ) -> pd.DataFrame:
        """
        Build a comprehensive triplicate-level view of all CT values.
        Returns DataFrame with one row per well, grouped by Sample+Target.
        """
        if data is None or data.empty:
            return pd.DataFrame()

        excluded_wells = excluded_wells or set()

        # Create base dataframe with all wells
        df = data[["Well", "Sample", "Target", "CT"]].copy()
        df["Excluded"] = df["Well"].isin(excluded_wells)

        # Calculate triplicate statistics per Sample+Target group
        group_stats = (
            df[~df["Excluded"]]
            .groupby(["Sample", "Target"])
            .agg({"CT": ["mean", "std", "count", "min", "max", list], "Well": list})
            .reset_index()
        )

        group_stats.columns = [
            "Sample",
            "Target",
            "Mean_CT",
            "SD",
            "n",
            "Min_CT",
            "Max_CT",
            "CT_Values",
            "Wells",
        ]

        # Calculate CV% and Range
        group_stats["CV_pct"] = np.where(
            group_stats["Mean_CT"] > 0,
            (group_stats["SD"] / group_stats["Mean_CT"]) * 100,
            0,
        )
        group_stats["Range"] = group_stats["Max_CT"] - group_stats["Min_CT"]

        # Determine health status for each group
        def get_health_status(row):
            issues = []
            severity = "ok"

            if row["n"] < 2:
                issues.append("Low n")
                severity = "warning"
            if row["CV_pct"] > QualityControl.CV_THRESHOLD * 100:
                issues.append(f"High CV ({row['CV_pct']:.1f}%)")
                severity = "warning"
            if row["Mean_CT"] > QualityControl.CT_HIGH_THRESHOLD:
                issues.append("High CT")
                severity = "warning"
            if row["Mean_CT"] < QualityControl.CT_LOW_THRESHOLD:
                issues.append("Low CT")
                severity = "warning"
            if row["Range"] > 1.0 and row["n"] >= 2:
                issues.append(f"High range ({row['Range']:.2f})")
                severity = "error" if row["Range"] > 2.0 else "warning"

            # Grubbs test for outliers (n >= 3)
            if row["n"] >= 3:
                ct_vals = np.array(row["CT_Values"])
                is_outlier, _ = QualityControl.grubbs_test(
                    ct_vals, QualityControl.GRUBBS_ALPHA
                )
                if is_outlier:
                    issues.append("Has outlier")
                    severity = "error"

            if not issues:
                return "OK", "ok"
            return "; ".join(issues), severity

        group_stats[["Status", "Severity"]] = group_stats.apply(
            lambda row: pd.Series(get_health_status(row)), axis=1
        )

        # Round numeric columns for display
        group_stats["Mean_CT"] = group_stats["Mean_CT"].round(2)
        group_stats["SD"] = group_stats["SD"].round(3)
        group_stats["CV_pct"] = group_stats["CV_pct"].round(1)
        group_stats["Range"] = group_stats["Range"].round(2)

        return group_stats

    @staticmethod
    def get_wells_for_triplicate(
        data: pd.DataFrame, sample: str, target: str
    ) -> pd.DataFrame:
        """
        Get all wells for a specific Sample+Target combination with detailed info.
        """
        if data is None or data.empty:
            return pd.DataFrame()

        wells = data[(data["Sample"] == sample) & (data["Target"] == target)].copy()

        if wells.empty:
            return pd.DataFrame()

        # Calculate stats for this group
        mean_ct = wells["CT"].mean()
        std_ct = wells["CT"].std() if len(wells) > 1 else 0

        # Mark outliers using Grubbs test
        wells["Is_Outlier"] = False
        if len(wells) >= 3:
            ct_vals = wells["CT"].values
            is_outlier, outlier_idx = QualityControl.grubbs_test(
                ct_vals, QualityControl.GRUBBS_ALPHA
            )
            if is_outlier and outlier_idx >= 0:
                wells.iloc[outlier_idx, wells.columns.get_loc("Is_Outlier")] = True

        # Calculate deviation from mean
        wells["Deviation"] = (wells["CT"] - mean_ct).round(3)
        wells["Z_Score"] = np.where(
            std_ct > 0, (wells["CT"] - mean_ct) / std_ct, 0
        ).round(2)

        # Determine individual well status
        def well_status(row):
            issues = []
            if row["CT"] > QualityControl.CT_HIGH_THRESHOLD:
                issues.append("High CT")
            if row["CT"] < QualityControl.CT_LOW_THRESHOLD:
                issues.append("Low CT")
            if row["Is_Outlier"]:
                issues.append("Grubbs outlier")
            if abs(row["Z_Score"]) > 2:
                issues.append(f"Z={row['Z_Score']:.1f}")
            return "; ".join(issues) if issues else "OK"

        wells["Well_Status"] = wells.apply(well_status, axis=1)

        return wells[
            [
                "Well",
                "Sample",
                "Target",
                "CT",
                "Deviation",
                "Z_Score",
                "Is_Outlier",
                "Well_Status",
            ]
        ]

    @staticmethod
    def get_qc_summary_stats(data: pd.DataFrame, excluded_wells: set = None) -> dict:
        """
        Calculate overall QC summary statistics for the dataset.
        """
        if data is None or data.empty:
            return {}

        excluded_wells = excluded_wells or set()
        active_data = data[~data["Well"].isin(excluded_wells)]

        # Basic counts
        total_wells = len(data)
        excluded_count = len(data[data["Well"].isin(excluded_wells)])
        active_wells = total_wells - excluded_count

        # CT distribution
        ct_mean = active_data["CT"].mean()
        ct_std = active_data["CT"].std()
        ct_min = active_data["CT"].min()
        ct_max = active_data["CT"].max()

        # Issue counts
        high_ct_count = len(
            active_data[active_data["CT"] > QualityControl.CT_HIGH_THRESHOLD]
        )
        low_ct_count = len(
            active_data[active_data["CT"] < QualityControl.CT_LOW_THRESHOLD]
        )

        # Triplicate health
        triplicate_stats = QualityControl.get_triplicate_data(data, excluded_wells)
        if not triplicate_stats.empty:
            total_triplicates = len(triplicate_stats)
            healthy_triplicates = len(
                triplicate_stats[triplicate_stats["Status"] == "OK"]
            )
            warning_triplicates = len(
                triplicate_stats[triplicate_stats["Severity"] == "warning"]
            )
            error_triplicates = len(
                triplicate_stats[triplicate_stats["Severity"] == "error"]
            )
            avg_cv = triplicate_stats["CV_pct"].mean()
            max_cv = triplicate_stats["CV_pct"].max()
        else:
            total_triplicates = healthy_triplicates = warning_triplicates = (
                error_triplicates
            ) = 0
            avg_cv = max_cv = 0

        # Unique counts
        n_samples = active_data["Sample"].nunique()
        n_genes = active_data["Target"].nunique()

        return {
            "total_wells": total_wells,
            "excluded_wells": excluded_count,
            "active_wells": active_wells,
            "ct_mean": round(ct_mean, 2) if pd.notna(ct_mean) else 0,
            "ct_std": round(ct_std, 2) if pd.notna(ct_std) else 0,
            "ct_min": round(ct_min, 2) if pd.notna(ct_min) else 0,
            "ct_max": round(ct_max, 2) if pd.notna(ct_max) else 0,
            "high_ct_count": high_ct_count,
            "low_ct_count": low_ct_count,
            "total_triplicates": total_triplicates,
            "healthy_triplicates": healthy_triplicates,
            "warning_triplicates": warning_triplicates,
            "error_triplicates": error_triplicates,
            "avg_cv_pct": round(avg_cv, 1) if pd.notna(avg_cv) else 0,
            "max_cv_pct": round(max_cv, 1) if pd.notna(max_cv) else 0,
            "n_samples": n_samples,
            "n_genes": n_genes,
            "health_score": round(healthy_triplicates / total_triplicates * 100, 1)
            if total_triplicates > 0
            else 0,
        }

    @staticmethod
    def suggest_exclusions(
        data: pd.DataFrame,
        sample: str,
        target: str,
        excluded_wells: set = None,
        strategy: str = "outlier",
    ) -> list:
        """
        Suggest wells to exclude based on different strategies.

        Strategies:
        - 'outlier': Exclude statistical outliers (Grubbs test)
        - 'worst': Exclude the well with highest deviation from mean
        - 'keep_best_2': Keep the 2 closest values, exclude others
        """
        excluded_wells = excluded_wells or set()
        wells_df = QualityControl.get_wells_for_triplicate(data, sample, target)

        if wells_df.empty:
            return []

        # Filter to active wells only
        active_wells = wells_df[~wells_df["Well"].isin(excluded_wells)]

        if len(active_wells) < 2:
            return []

        suggestions = []

        if strategy == "outlier":
            # Suggest Grubbs outliers
            outliers = active_wells[active_wells["Is_Outlier"]]
            suggestions = outliers["Well"].tolist()

        elif strategy == "worst":
            # Suggest well with highest absolute deviation
            if len(active_wells) > 2:
                worst_idx = active_wells["Deviation"].abs().idxmax()
                suggestions = [active_wells.loc[worst_idx, "Well"]]

        elif strategy == "keep_best_2":
            # Keep the 2 closest values to the median
            if len(active_wells) > 2:
                median_ct = active_wells["CT"].median()
                active_wells_sorted = active_wells.copy()
                active_wells_sorted["Dist_to_Median"] = abs(
                    active_wells_sorted["CT"] - median_ct
                )
                active_wells_sorted = active_wells_sorted.sort_values("Dist_to_Median")
                # Exclude all but the best 2
                to_exclude = active_wells_sorted.iloc[2:]["Well"].tolist()
                suggestions = to_exclude

        return suggestions

    @staticmethod
    def grubbs_test(values: np.ndarray, alpha: float = 0.05) -> Tuple[bool, int]:
        n = len(values)
        if n < 3:
            return False, -1

        mean_val = np.mean(values)
        std_val = np.std(values, ddof=1)

        if std_val == 0:
            return False, -1

        g_scores = np.abs(values - mean_val) / std_val
        max_idx = np.argmax(g_scores)
        g_stat = g_scores[max_idx]

        t_crit = stats.t.ppf(1 - alpha / (2 * n), n - 2)
        g_crit = ((n - 1) / np.sqrt(n)) * np.sqrt(t_crit**2 / (n - 2 + t_crit**2))

        return g_stat > g_crit, int(max_idx)

    @staticmethod
    def detect_outliers(data: pd.DataFrame, hk_gene: str = None) -> pd.DataFrame:
        if data is None or data.empty:
            return pd.DataFrame()

        qc_df = data[["Well", "Sample", "Target", "CT"]].copy()

        ct_high = qc_df["CT"] > QualityControl.CT_HIGH_THRESHOLD
        ct_low = qc_df["CT"] < QualityControl.CT_LOW_THRESHOLD

        high_ct_issue = f"CT > {QualityControl.CT_HIGH_THRESHOLD} (low expression)"
        low_ct_issue = f"CT < {QualityControl.CT_LOW_THRESHOLD} (unusually high)"

        qc_df["Issues"] = pd.DataFrame(
            {
                "high": ct_high.map(lambda x: high_ct_issue if x else ""),
                "low": ct_low.map(lambda x: low_ct_issue if x else ""),
            }
        ).apply(lambda row: "; ".join([x for x in row if x]) or "OK", axis=1)

        qc_df["Severity"] = np.where(ct_high | ct_low, "warning", "ok")
        qc_df["Flagged"] = ct_high | ct_low

        cv_stats = (
            data.groupby(["Sample", "Target"])["CT"]
            .agg(["mean", "std", "count"])
            .reset_index()
        )
        cv_stats["cv"] = np.where(
            (cv_stats["mean"] > 0) & (cv_stats["count"] > 1),
            cv_stats["std"] / cv_stats["mean"],
            0,
        )
        high_cv_groups = cv_stats[cv_stats["cv"] > QualityControl.CV_THRESHOLD][
            ["Sample", "Target", "cv"]
        ]

        if not high_cv_groups.empty:
            qc_df = qc_df.merge(high_cv_groups, on=["Sample", "Target"], how="left")
            has_high_cv = qc_df["cv"].notna()

            cv_issue = qc_df["cv"].apply(
                lambda x: f"CV={x:.1%} (high variability)" if pd.notna(x) else ""
            )
            qc_df.loc[has_high_cv, "Issues"] = qc_df.loc[has_high_cv].apply(
                lambda r: f"{r['Issues']}; {cv_issue[r.name]}"
                if r["Issues"] != "OK"
                else cv_issue[r.name],
                axis=1,
            )
            qc_df.loc[has_high_cv, "Severity"] = "warning"
            qc_df.loc[has_high_cv, "Flagged"] = True
            qc_df = qc_df.drop(columns=["cv"])

        grubbs_outliers = set()
        for (sample, target), group in data.groupby(["Sample", "Target"]):
            if len(group) >= 3:
                ct_vals = group["CT"].values
                is_outlier, outlier_idx = QualityControl.grubbs_test(
                    ct_vals, QualityControl.GRUBBS_ALPHA
                )
                if is_outlier:
                    outlier_well = group.iloc[outlier_idx]["Well"]
                    grubbs_outliers.add(outlier_well)

        if grubbs_outliers:
            grubbs_mask = qc_df["Well"].isin(grubbs_outliers)

            def add_grubbs_issue(current_issue):
                grubbs_issue = "Grubbs outlier"
                return (
                    f"{current_issue}; {grubbs_issue}"
                    if current_issue != "OK"
                    else grubbs_issue
                )

            qc_df.loc[grubbs_mask, "Issues"] = qc_df.loc[grubbs_mask, "Issues"].apply(
                add_grubbs_issue
            )
            qc_df.loc[grubbs_mask, "Severity"] = "error"
            qc_df.loc[grubbs_mask, "Flagged"] = True

        if hk_gene:
            hk_data = data[data["Target"] == hk_gene]
            if not hk_data.empty:
                hk_by_sample = hk_data.groupby("Sample")["CT"].mean()
                overall_hk_mean = hk_by_sample.mean()

                deviations = (hk_by_sample - overall_hk_mean).abs()
                flagged_samples = deviations[
                    deviations > QualityControl.HK_VARIATION_THRESHOLD
                ]

                if not flagged_samples.empty:
                    deviation_map = flagged_samples.to_dict()
                    hk_mask = (qc_df["Target"] == hk_gene) & (
                        qc_df["Sample"].isin(flagged_samples.index)
                    )

                    def add_hk_issue(row):
                        dev = deviation_map.get(row["Sample"], 0)
                        hk_issue = f"HK deviation={dev:.2f}"
                        return (
                            f"{row['Issues']}; {hk_issue}"
                            if row["Issues"] != "OK"
                            else hk_issue
                        )

                    qc_df.loc[hk_mask, "Issues"] = qc_df.loc[hk_mask].apply(
                        add_hk_issue, axis=1
                    )
                    qc_df.loc[hk_mask, "Severity"] = "error"
                    qc_df.loc[hk_mask, "Flagged"] = True

        return qc_df

    @staticmethod
    def get_replicate_stats(data: pd.DataFrame) -> pd.DataFrame:
        if data is None or data.empty:
            return pd.DataFrame()

        stats = (
            data.groupby(["Sample", "Target"])["CT"]
            .agg(["mean", "std", "count"])
            .reset_index()
        )
        stats.columns = ["Sample", "Target", "Mean CT", "SD", "n"]

        stats["SD"] = stats["SD"].fillna(0)
        stats["CV%"] = np.where(
            stats["Mean CT"] > 0, (stats["SD"] / stats["Mean CT"]) * 100, 0
        )

        stats["Status"] = np.select(
            [stats["Mean CT"] < 10, stats["Mean CT"] > 35, stats["CV%"] > 5],
            ["Check Signal", "Low Expression", "High CV"],
            default="OK",
        )

        stats["Mean CT"] = stats["Mean CT"].round(2)
        stats["SD"] = stats["SD"].round(3)
        stats["CV%"] = stats["CV%"].round(1)
        stats["n"] = stats["n"].astype(int)

        return stats[["Sample", "Target", "n", "Mean CT", "SD", "CV%", "Status"]]

    @staticmethod
    def create_plate_heatmap(
        data: pd.DataFrame, value_col: str = "CT", excluded_wells: set = None
    ) -> go.Figure:
        if data is None or data.empty:
            return go.Figure()

        excluded_wells = excluded_wells or set()

        rows = list("ABCDEFGH")
        cols = list(range(1, 13))

        plate_values = np.full((8, 12), np.nan)
        plate_text = [["" for _ in range(12)] for _ in range(8)]
        plate_colors = [[0 for _ in range(12)] for _ in range(8)]

        for _, row in data.iterrows():
            well = row["Well"]
            if len(well) >= 2:
                well_row = well[0].upper()
                try:
                    well_col = int(well[1:])
                except ValueError:
                    continue

                if well_row in rows and 1 <= well_col <= 12:
                    r_idx = rows.index(well_row)
                    c_idx = well_col - 1

                    ct_val = row[value_col]
                    plate_values[r_idx, c_idx] = ct_val

                    sample_short = str(row["Sample"])[:10]
                    target_short = str(row["Target"])[:8]
                    excluded_marker = " [X]" if well in excluded_wells else ""
                    plate_text[r_idx][c_idx] = (
                        f"{well}{excluded_marker}<br>{sample_short}<br>{target_short}<br>CT: {ct_val:.1f}"
                    )

                    if well in excluded_wells:
                        plate_colors[r_idx][c_idx] = -1

        fig = go.Figure(
            data=go.Heatmap(
                z=plate_values,
                x=[str(c) for c in cols],
                y=rows,
                text=plate_text,
                hoverinfo="text",
                colorscale=[[0, "#2ecc71"], [0.5, "#f1c40f"], [1, "#e74c3c"]],
                zmin=15,
                zmax=40,
                colorbar=dict(title="CT Value"),
            )
        )

        for r_idx, row_letter in enumerate(rows):
            for c_idx, col_num in enumerate(cols):
                well_name = f"{row_letter}{col_num}"
                if well_name in excluded_wells:
                    fig.add_annotation(
                        x=str(col_num),
                        y=row_letter,
                        text="X",
                        showarrow=False,
                        font=dict(size=20, color="red"),
                        opacity=0.8,
                    )

        fig.update_layout(
            title="96-Well Plate Overview",
            xaxis=dict(title="Column", side="top", dtick=1),
            yaxis=dict(title="Row", autorange="reversed"),
            height=400,
            width=800,
        )

        return fig


# ==================== ANALYSIS ENGINE ====================
class AnalysisEngine:
    @staticmethod
    def calculate_ddct(
        data: pd.DataFrame,
        hk_gene: str,
        ref_sample: str,
        excluded_wells: set,
        excluded_samples: set,
        sample_mapping: dict,
    ) -> pd.DataFrame:
        """Gene-by-gene ΔΔCt calculation with housekeeping normalization"""

        # Filter data
        data = data[
            ~data["Well"].isin(excluded_wells) & ~data["Sample"].isin(excluded_samples)
        ].copy()

        # Apply sample name mapping
        data["Condition"] = data["Sample"].map(
            lambda x: sample_mapping.get(x, {}).get("condition", x)
        )
        data["Group"] = data["Sample"].map(
            lambda x: sample_mapping.get(x, {}).get("group", "Treatment")
        )

        results = []

        # Process each target gene separately (exclude housekeeping)
        for target in data["Target"].unique():
            if target.upper() in [hk_gene.upper(), "ACTIN", "B-ACTIN", "GAPDH", "ACTB"]:
                continue

            target_data = data[data["Target"] == target]

            for condition in target_data["Condition"].unique():
                cond_data = target_data[target_data["Condition"] == condition]

                if len(cond_data) == 0:
                    continue

                hk_data = data[
                    (data["Condition"] == condition) & (data["Target"] == hk_gene)
                ]
                if len(hk_data) == 0:
                    continue

                target_ct_values = cond_data["CT"].values
                hk_ct_values = hk_data["CT"].values

                if len(target_ct_values) == 0 or len(hk_ct_values) == 0:
                    continue

                target_ct_mean = target_ct_values.mean()
                hk_ct_mean = hk_ct_values.mean()
                delta_ct = target_ct_mean - hk_ct_mean

                # Get reference ΔCt (ref_sample)
                ref_target = target_data[target_data["Condition"] == ref_sample]
                ref_hk = data[
                    (data["Condition"] == ref_sample) & (data["Target"] == hk_gene)
                ]

                if len(ref_target) > 0 and len(ref_hk) > 0:
                    ref_delta_ct = ref_target["CT"].mean() - ref_hk["CT"].mean()
                else:
                    ref_delta_ct = 0

                ddct = delta_ct - ref_delta_ct
                rel_expr = 2 ** (-ddct)

                target_sd = target_ct_values.std() if len(target_ct_values) > 1 else 0
                hk_sd = hk_ct_values.std() if len(hk_ct_values) > 1 else 0
                n_target = len(target_ct_values)
                n_hk = len(hk_ct_values)

                target_sem = target_sd / np.sqrt(n_target) if n_target > 1 else 0
                hk_sem = hk_sd / np.sqrt(n_hk) if n_hk > 1 else 0
                combined_sem = np.sqrt(target_sem**2 + hk_sem**2)

                sem = combined_sem * rel_expr * np.log(2)

                # Get original sample name and group
                original_sample = cond_data["Sample"].iloc[0]
                group = sample_mapping.get(original_sample, {}).get(
                    "group", "Treatment"
                )

                results.append(
                    {
                        "Target": target,
                        "Condition": condition,
                        "Original_Sample": original_sample,
                        "Group": group,
                        "n_replicates": n_target,
                        "n_hk_replicates": n_hk,
                        "Target_Ct_Mean": target_ct_mean,
                        "Target_Ct_SD": target_sd,
                        "HK_Ct_Mean": hk_ct_mean,
                        "HK_Ct_SD": hk_sd,
                        "Delta_Ct": delta_ct,
                        "Delta_Delta_Ct": ddct,
                        "Relative_Expression": rel_expr,
                        "SEM": sem,
                        "Fold_Change": rel_expr,
                    }
                )

        return pd.DataFrame(results)

    @staticmethod
    def calculate_statistics(
        processed: pd.DataFrame,
        compare_condition: str,
        compare_condition_2: str = None,
        raw_data: pd.DataFrame = None,
        hk_gene: str = None,
        sample_mapping: dict = None,
    ) -> pd.DataFrame:
        """Two-tailed Welch's t-test comparing each condition to compare_condition (and optionally compare_condition_2)"""

        # Use session_state fallbacks
        raw_data = raw_data if raw_data is not None else st.session_state.get("data")
        hk_gene = hk_gene if hk_gene is not None else st.session_state.get("hk_gene")
        sample_mapping = (
            sample_mapping
            if sample_mapping is not None
            else st.session_state.get("sample_mapping", {})
        )

        if raw_data is None or hk_gene is None:
            return processed

        results = processed.copy()
        results["p_value"] = np.nan
        results["significance"] = ""

        # Add second p-value columns if compare_condition_2 is provided
        if compare_condition_2:
            results["p_value_2"] = np.nan
            results["significance_2"] = ""

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

            # FIRST COMPARISON: compare_condition
            ref_vals = rel_expr.get(compare_condition, np.array([]))
            if ref_vals.size >= 1:
                for cond, vals in rel_expr.items():
                    if cond == compare_condition or vals.size == 0:
                        continue

                    try:
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore", RuntimeWarning)
                            if ref_vals.size >= 2 and vals.size >= 2:
                                _, p_val = stats.ttest_ind(
                                    ref_vals, vals, equal_var=False
                                )
                            elif vals.size == 1 and ref_vals.size >= 2:
                                _, p_val = stats.ttest_1samp(ref_vals, vals[0])
                            elif ref_vals.size == 1 and vals.size >= 2:
                                _, p_val = stats.ttest_1samp(vals, ref_vals[0])
                            else:
                                p_val = np.nan
                    except (ValueError, TypeError) as e:
                        p_val = np.nan

                    mask = (results["Target"] == target) & (
                        results["Condition"] == cond
                    )
                    results.loc[mask, "p_value"] = p_val

                    if not np.isnan(p_val):
                        if p_val < 0.001:
                            results.loc[mask, "significance"] = "***"
                        elif p_val < 0.01:
                            results.loc[mask, "significance"] = "**"
                        elif p_val < 0.05:
                            results.loc[mask, "significance"] = "*"

            # SECOND COMPARISON: compare_condition_2 (if provided)
            if compare_condition_2:
                ref_vals_2 = rel_expr.get(compare_condition_2, np.array([]))
                if ref_vals_2.size >= 1:
                    for cond, vals in rel_expr.items():
                        if cond == compare_condition_2 or vals.size == 0:
                            continue

                        try:
                            with warnings.catch_warnings():
                                warnings.simplefilter("ignore", RuntimeWarning)
                                if ref_vals_2.size >= 2 and vals.size >= 2:
                                    _, p_val_2 = stats.ttest_ind(
                                        ref_vals_2, vals, equal_var=False
                                    )
                                elif vals.size == 1 and ref_vals_2.size >= 2:
                                    _, p_val_2 = stats.ttest_1samp(ref_vals_2, vals[0])
                                elif ref_vals_2.size == 1 and vals.size >= 2:
                                    _, p_val_2 = stats.ttest_1samp(vals, ref_vals_2[0])
                                else:
                                    p_val_2 = np.nan
                        except (ValueError, TypeError) as e:
                            p_val_2 = np.nan

                        mask = (results["Target"] == target) & (
                            results["Condition"] == cond
                        )
                        results.loc[mask, "p_value_2"] = p_val_2

                        if not np.isnan(p_val_2):
                            if p_val_2 < 0.001:
                                results.loc[mask, "significance_2"] = "###"
                            elif p_val_2 < 0.01:
                                results.loc[mask, "significance_2"] = "##"
                            elif p_val_2 < 0.05:
                                results.loc[mask, "significance_2"] = "#"

        results = AnalysisEngine._apply_fdr_correction(
            results, "p_value", "p_value_fdr", "significance_fdr", "*"
        )
        if compare_condition_2:
            results = AnalysisEngine._apply_fdr_correction(
                results, "p_value_2", "p_value_fdr_2", "significance_fdr_2", "#"
            )

        return results

    @staticmethod
    def _apply_fdr_correction(
        results: pd.DataFrame, p_col: str, fdr_col: str, sig_col: str, marker: str
    ) -> pd.DataFrame:
        valid_mask = results[p_col].notna()
        valid_pvals = results.loc[valid_mask, p_col].values

        if len(valid_pvals) == 0:
            results[fdr_col] = np.nan
            results[sig_col] = ""
            return results

        n = len(valid_pvals)
        sorted_idx = np.argsort(valid_pvals)
        sorted_pvals = valid_pvals[sorted_idx]

        fdr_vals = np.zeros(n)
        for i in range(n):
            rank = i + 1
            fdr_vals[i] = sorted_pvals[i] * n / rank

        for i in range(n - 2, -1, -1):
            fdr_vals[i] = min(fdr_vals[i], fdr_vals[i + 1])

        fdr_vals = np.clip(fdr_vals, 0, 1)

        unsorted_fdr = np.zeros(n)
        unsorted_fdr[sorted_idx] = fdr_vals

        results[fdr_col] = np.nan
        results.loc[valid_mask, fdr_col] = unsorted_fdr

        results[sig_col] = ""
        for i, (idx, fdr) in enumerate(zip(results.index[valid_mask], unsorted_fdr)):
            if fdr < 0.001:
                results.loc[idx, sig_col] = marker * 3
            elif fdr < 0.01:
                results.loc[idx, sig_col] = marker * 2
            elif fdr < 0.05:
                results.loc[idx, sig_col] = marker

        return results

    @staticmethod
    def run_full_analysis(
        ref_sample_key: str, compare_sample_key: str, compare_sample_key_2: str = None
    ):
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

            ref_condition = mapping.get(ref_sample_key, {}).get(
                "condition", ref_sample_key
            )
            cmp_condition = mapping.get(compare_sample_key, {}).get(
                "condition", compare_sample_key
            )
            cmp_condition_2 = None
            if compare_sample_key_2:
                cmp_condition_2 = mapping.get(compare_sample_key_2, {}).get(
                    "condition", compare_sample_key_2
                )

            st.session_state.analysis_ref_condition = ref_condition
            st.session_state.analysis_cmp_condition = cmp_condition
            st.session_state.analysis_cmp_condition_2 = cmp_condition_2

            msg = f"Running full analysis using reference '{ref_condition}' and comparison '{cmp_condition}'"
            if cmp_condition_2:
                msg += f" + secondary comparison '{cmp_condition_2}'"

            with st.spinner(msg + "..."):
                # --- ΔΔCt calculation ---
                processed_df = AnalysisEngine.calculate_ddct(
                    data,
                    hk_gene,
                    ref_condition,
                    st.session_state.get("excluded_wells", set()),
                    st.session_state.get("excluded_samples", set()),
                    mapping,
                )

                if processed_df is None or processed_df.empty:
                    st.warning(
                        "⚠️ No ΔΔCt results produced. Check mapping and housekeeping gene."
                    )
                    return False

                # --- Statistical test ---
                try:
                    processed_with_stats = AnalysisEngine.calculate_statistics(
                        processed_df,
                        cmp_condition,
                        cmp_condition_2,
                        raw_data=data,
                        hk_gene=hk_gene,
                        sample_mapping=mapping,
                    )
                except TypeError:
                    processed_with_stats = AnalysisEngine.calculate_statistics(
                        processed_df, cmp_condition
                    )

                # --- Organize data for graphs ---
                gene_dict = {}
                if "Target" in processed_with_stats.columns:
                    for gene in processed_with_stats["Target"].unique():
                        gene_df = processed_with_stats[
                            processed_with_stats["Target"] == gene
                        ].copy()
                        gene_dict[gene] = gene_df.reset_index(drop=True)
                else:
                    gene_dict = {"results": processed_with_stats.reset_index(drop=True)}

                st.session_state.processed_data = gene_dict

            st.success(
                "✅ Full analysis complete. Go to the Graphs tab to visualize results."
            )
            return True

        except Exception as e:
            st.error(f"❌ Analysis failed: {e}")
            return False


# ==================== GRAPH GENERATOR ====================
import textwrap


class GraphGenerator:
    @staticmethod
    def _wrap_text(text: str, width: int = 15) -> str:
        """Wrap text for x-axis labels"""
        wrapped = textwrap.fill(text, width=width)
        return wrapped

    @staticmethod
    def create_gene_graph(
        data: pd.DataFrame,
        gene: str,
        settings: dict,
        efficacy_config: dict = None,
        sample_order: list = None,
        per_sample_overrides: dict = None,
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
        if "Relative_Expression" not in gene_data.columns:
            if "Fold_Change" in gene_data.columns:
                gene_data["Relative_Expression"] = gene_data["Fold_Change"]
            else:
                st.error(
                    f"Missing Relative_Expression or Fold_Change column for {gene}"
                )
                fig = go.Figure()
                fig.add_annotation(
                    text=f"Missing data columns for {gene}", showarrow=False
                )
                return fig

        if "SEM" not in gene_data.columns:
            gene_data["SEM"] = 0

        # FIXED: Use sample_order from mapping and deduplicate conditions while preserving order
        if sample_order:
            # Convert sample names to condition names using mapping
            mapping = st.session_state.get("sample_mapping", {})
            condition_order = []
            seen_conditions = set()  # Track which conditions we've already added

            for sample in sample_order:
                # Only include samples that are marked as 'include'
                if mapping.get(sample, {}).get("include", True):
                    cond = mapping.get(sample, {}).get("condition", sample)
                    # Only add if this condition exists in the data AND we haven't added it yet
                    if (
                        cond in gene_data["Condition"].unique()
                        and cond not in seen_conditions
                    ):
                        condition_order.append(cond)
                        seen_conditions.add(cond)

            # Add any conditions not in order (shouldn't happen, but safety)
            for cond in gene_data["Condition"].unique():
                if cond not in seen_conditions:
                    condition_order.append(cond)
                    seen_conditions.add(cond)

            # Apply categorical ordering
            gene_data["Condition"] = pd.Categorical(
                gene_data["Condition"], categories=condition_order, ordered=True
            )
            gene_data = gene_data.sort_values("Condition")
        else:
            # Fallback: sort by appearance order
            gene_data = gene_data.sort_values("Condition")

        # Reset index to ensure proper sequential indexing
        gene_data_indexed = gene_data.reset_index(drop=True)

        # Store condition names for labels
        condition_names = gene_data_indexed["Condition"].tolist()
        n_bars = len(gene_data_indexed)

        # Get colors - White/Medium Grey for controls, Grey base for treatments
        bar_colors = []

        for idx, row in gene_data_indexed.iterrows():
            condition = row["Condition"]
            group = row.get("Group", "Treatment")

            custom_key = f"{gene}_{condition}"
            if custom_key in settings.get("bar_colors_per_sample", {}):
                bar_colors.append(settings["bar_colors_per_sample"][custom_key])
            elif group in DEFAULT_GROUP_COLORS:
                bar_colors.append(DEFAULT_GROUP_COLORS[group])
            else:
                default_color = settings.get("bar_colors", {}).get(gene, "#D3D3D3")
                bar_colors.append(default_color)

        # Create figure
        fig = go.Figure()

        # Error bars
        error_array = (
            gene_data_indexed["SEM"] * settings.get("error_multiplier", 1.96)
        ).values

        # Per-bar settings for individual control
        gene_bar_settings = st.session_state.get(f"{gene}_bar_settings", {})

        # Check global show/hide for this gene
        show_error_global = settings.get("show_error", True)
        show_sig_global = settings.get("show_significance", True)

        # Build error visibility array
        error_visible_array = []

        for idx in range(n_bars):
            row = gene_data_indexed.iloc[idx]
            condition = row["Condition"]
            bar_key = f"{gene}_{condition}"

            # Get individual bar settings (default to True)
            bar_config = gene_bar_settings.get(
                bar_key, {"show_sig": True, "show_err": True}
            )

            # ERROR BARS: Both global AND individual must be True
            if show_error_global and bar_config.get("show_err", True):
                error_visible_array.append(error_array[idx])
            else:
                error_visible_array.append(0)

        # Add bar trace with UPPER-ONLY error bars
        # CRITICAL: Use numeric x-values (indices) for proper positioning
        fig.add_trace(
            go.Bar(
                x=list(range(n_bars)),  # Use numeric indices 0, 1, 2, ... n_bars-1
                y=gene_data_indexed["Relative_Expression"],
                error_y=dict(
                    type="data",
                    array=error_visible_array,
                    arrayminus=[0] * n_bars,  # NO LOWER ERROR BARS
                    visible=True,
                    thickness=2,
                    width=4,
                    color="rgba(0,0,0,0.5)",
                    symmetric=False,
                ),
                marker=dict(
                    color=bar_colors,
                    line=dict(
                        width=settings.get("marker_line_width", 1), color="black"
                    ),
                    opacity=settings.get("bar_opacity", 0.95),
                ),
                showlegend=False,
            )
        )

        # Calculate y-axis range FIRST (needed for absolute positioning of significance)
        max_y_value = gene_data_indexed["Relative_Expression"].max()
        max_error = error_array.max() if len(error_array) > 0 else 0
        y_max_auto = (
            max_y_value + max_error + (max_y_value * 0.15)
        )  # Add 15% padding for stars

        # FIXED ABSOLUTE SPACING for dual symbols (in data units)
        # This ensures consistent spacing across all bars regardless of height
        fixed_symbol_spacing = y_max_auto * 0.05  # 5% of y-axis as fixed spacing unit

        # Add significance symbols - aligned with bars (DUAL SUPPORT with absolute positioning)
        for idx in range(n_bars):
            row = gene_data_indexed.iloc[idx]
            condition = row["Condition"]
            bar_key = f"{gene}_{condition}"
            bar_config = gene_bar_settings.get(
                bar_key, {"show_sig": True, "show_err": True}
            )

            # Get both significance values
            sig_1 = row.get("significance", "")
            sig_2 = row.get("significance_2", "")

            # Calculate base y position (top of error bar)
            bar_height = row["Relative_Expression"]
            error_bar_height = error_visible_array[idx]
            base_y_position = bar_height + error_bar_height

            # Font sizes - asterisk at normal size, hashtag reduced to match visually
            asterisk_font_size = 16
            hashtag_font_size = 10  # Reduced by 6 to match asterisk visual size

            # Check if we need to show significance
            if show_sig_global and bar_config.get("show_sig", True):
                symbols_to_show = []
                font_sizes = []

                # Add first significance (asterisks)
                if sig_1 in ["*", "**", "***"]:
                    symbols_to_show.append(sig_1)
                    font_sizes.append(asterisk_font_size)

                # Add second significance (hashtags)
                if sig_2 in ["#", "##", "###"]:
                    symbols_to_show.append(sig_2)
                    font_sizes.append(hashtag_font_size)

                # Display symbols with FIXED ABSOLUTE SPACING
                if len(symbols_to_show) == 2:
                    # Two symbols: stack with FIXED absolute spacing
                    # Bottom symbol (asterisk) - positioned above error bar
                    fig.add_annotation(
                        x=idx,
                        y=base_y_position + (fixed_symbol_spacing * 0.2),
                        text=symbols_to_show[0],
                        showarrow=False,
                        font=dict(size=font_sizes[0], color="black", family="Arial"),
                        xref="x",
                        yref="y",
                        xanchor="center",
                        yanchor="bottom",
                    )

                    # Top symbol (hashtag) - FIXED absolute distance above bottom symbol
                    # The vertical gap is always fixed_symbol_spacing regardless of bar height
                    fig.add_annotation(
                        x=idx,
                        y=base_y_position
                        + (fixed_symbol_spacing * 0.2)
                        + fixed_symbol_spacing,
                        text=symbols_to_show[1],
                        showarrow=False,
                        font=dict(size=font_sizes[1], color="black", family="Arial"),
                        xref="x",
                        yref="y",
                        xanchor="center",
                        yanchor="bottom",
                    )

                elif len(symbols_to_show) == 1:
                    # Single symbol - positioned just above error bar
                    fig.add_annotation(
                        x=idx,
                        y=base_y_position + (fixed_symbol_spacing * 0.2),
                        text=symbols_to_show[0],
                        showarrow=False,
                        font=dict(size=font_sizes[0], color="black", family="Arial"),
                        xref="x",
                        yref="y",
                        xanchor="center",
                        yanchor="bottom",
                    )

        # Custom y-axis label with bold red gene name
        y_label_html = f"Relative <b style='color:red;'>{gene}</b> Expression Level"

        # Y-axis configuration (y_max_auto already calculated above)
        y_axis_config = dict(
            title=dict(
                text=y_label_html,
                font=dict(size=settings.get(f"{gene}_ylabel_size", 14)),
            ),
            showgrid=False,
            zeroline=True,
            zerolinewidth=1,
            zerolinecolor="black",
            range=[0, y_max_auto],
            fixedrange=False,
        )

        # ---
        if settings.get("y_log_scale"):
            y_axis_config["type"] = "log"
            y_axis_config.pop("range", None)

        # Manual range override if user specified
        if settings.get("y_min") is not None or settings.get("y_max") is not None:
            y_range = []
            y_range.append(settings.get("y_min", 0))
            y_range.append(settings.get("y_max", y_max_auto))
            y_axis_config["range"] = y_range

        # Get gene-specific settings
        gene_bar_gap = settings.get(f"{gene}_bar_gap", settings.get("bar_gap", 0.15))
        gene_margins = settings.get(
            f"{gene}_margins", {"l": 80, "r": 80, "t": 100, "b": 100}
        )
        gene_bg_color = settings.get(
            f"{gene}_bg_color", settings.get("plot_bgcolor", "#FFFFFF")
        )
        gene_tick_size = settings.get(f"{gene}_tick_size", 12)

        # Wrap x-axis labels - all of them
        wrapped_labels = [
            GraphGenerator._wrap_text(str(cond), 15) for cond in condition_names
        ]

        # P-VALUE LEGEND - Support dual comparison
        legend_text = "<b>Significance:</b>  * p<0.05  ** p<0.01  *** p<0.001"

        # Check if there's a second p-value comparison in the data
        if (
            "significance_2" in gene_data_indexed.columns
            and gene_data_indexed["significance_2"].notna().any()
        ):
            legend_text += (
                "<br><b>2nd Comparison:</b>  # p<0.05  ## p<0.01  ### p<0.001"
            )

        fig.add_annotation(
            text=legend_text,
            xref="paper",
            yref="paper",
            x=1.0,
            y=-0.15,
            xanchor="right",
            yanchor="top",
            showarrow=False,
            font=dict(size=12, color="#666666", family="Arial"),
            bgcolor="rgba(255,255,255,0.9)",
            bordercolor="#CCCCCC",
            borderwidth=1,
            borderpad=4,
        )

        fig.update_layout(
            title=dict(
                text=f"{gene} Expression",
                font=dict(
                    size=settings.get("title_size", 20), family="Arial", color="#333333"
                ),
                x=0.5,
                xanchor="center",
                y=0.98,
                yanchor="top",
            ),
            xaxis=dict(
                title=None,
                showgrid=False,
                zeroline=False,
                tickmode="array",
                tickvals=list(range(n_bars)),
                ticktext=wrapped_labels,
                tickfont=dict(size=gene_tick_size),
                tickangle=0,
                showline=False,
                mirror=False,
                side="bottom",
                range=[-0.5, n_bars - 0.5],
            ),
            yaxis=y_axis_config,
            template=settings.get("color_scheme", "plotly_white"),
            font=dict(size=settings.get("font_size", 14)),
            height=settings.get("figure_height", 600),
            width=settings.get("figure_width", 1000),
            bargap=gene_bar_gap,
            showlegend=settings.get("show_legend", False),
            plot_bgcolor=gene_bg_color,
            paper_bgcolor="#FFFFFF",
            margin=dict(
                l=gene_margins.get("l", 80),
                r=gene_margins.get("r", 80),
                t=gene_margins.get("t", 100),
                b=gene_margins.get("b", 120),
            ),
        )

        return fig


# ==================== REPORT GENERATOR (PPT) ====================
class ReportGenerator:
    SLIDE_WIDTH_INCHES = 13.333
    SLIDE_HEIGHT_INCHES = 7.5

    @staticmethod
    def create_presentation(
        graphs: Dict[str, go.Figure],
        processed_data: Dict[str, pd.DataFrame],
        analysis_params: dict,
        layout: str = "one_per_slide",
        include_title_slide: bool = True,
        include_summary: bool = True,
    ) -> bytes:
        try:
            from pptx import Presentation
            from pptx.util import Inches, Pt
            from pptx.dml.color import RgbColor
            from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
        except ImportError:
            raise ImportError(
                "python-pptx is required. Install with: pip install python-pptx"
            )

        prs = Presentation()
        prs.slide_width = Inches(ReportGenerator.SLIDE_WIDTH_INCHES)
        prs.slide_height = Inches(ReportGenerator.SLIDE_HEIGHT_INCHES)

        blank_layout = prs.slide_layouts[6]

        if include_title_slide:
            ReportGenerator._add_title_slide(prs, analysis_params)

        if layout == "one_per_slide":
            for gene, fig in graphs.items():
                gene_data = processed_data.get(gene)
                ReportGenerator._add_gene_slide(prs, blank_layout, gene, fig, gene_data)
        elif layout == "two_per_slide":
            gene_list = list(graphs.items())
            for i in range(0, len(gene_list), 2):
                pair = gene_list[i : i + 2]
                ReportGenerator._add_dual_gene_slide(
                    prs, blank_layout, pair, processed_data
                )
        elif layout == "grid":
            ReportGenerator._add_grid_slide(prs, blank_layout, graphs)

        if include_summary:
            ReportGenerator._add_summary_slide(
                prs, blank_layout, processed_data, analysis_params
            )

        output = io.BytesIO()
        prs.save(output)
        output.seek(0)
        return output.getvalue()

    @staticmethod
    def _add_title_slide(prs, analysis_params: dict):
        from pptx.util import Inches, Pt
        from pptx.enum.text import PP_ALIGN

        slide = prs.slides.add_slide(prs.slide_layouts[6])

        title_box = slide.shapes.add_textbox(
            Inches(0.5), Inches(2.5), Inches(12.33), Inches(1)
        )
        title_frame = title_box.text_frame
        title_para = title_frame.paragraphs[0]
        title_para.text = "qPCR Gene Expression Analysis"
        title_para.font.size = Pt(44)
        title_para.font.bold = True
        title_para.alignment = PP_ALIGN.CENTER

        subtitle_box = slide.shapes.add_textbox(
            Inches(0.5), Inches(3.8), Inches(12.33), Inches(0.5)
        )
        subtitle_frame = subtitle_box.text_frame
        subtitle_para = subtitle_frame.paragraphs[0]
        efficacy = analysis_params.get("Efficacy_Type", "Analysis")
        subtitle_para.text = f"{efficacy} Efficacy Study"
        subtitle_para.font.size = Pt(28)
        subtitle_para.alignment = PP_ALIGN.CENTER

        info_box = slide.shapes.add_textbox(
            Inches(0.5), Inches(5.5), Inches(12.33), Inches(1)
        )
        info_frame = info_box.text_frame
        info_para = info_frame.paragraphs[0]
        date = analysis_params.get("Date", datetime.now().strftime("%Y-%m-%d"))
        hk = analysis_params.get("Housekeeping_Gene", "N/A")
        ref = analysis_params.get("Reference_Sample", "N/A")
        info_para.text = f"Date: {date}  |  HK Gene: {hk}  |  Reference: {ref}"
        info_para.font.size = Pt(14)
        info_para.alignment = PP_ALIGN.CENTER

    @staticmethod
    def _add_gene_slide(
        prs, layout, gene: str, fig: go.Figure, gene_data: pd.DataFrame = None
    ):
        from pptx.util import Inches, Pt
        from pptx.enum.text import PP_ALIGN

        slide = prs.slides.add_slide(layout)

        title_box = slide.shapes.add_textbox(
            Inches(0.3), Inches(0.2), Inches(12.73), Inches(0.6)
        )
        title_frame = title_box.text_frame
        title_para = title_frame.paragraphs[0]
        title_para.text = f"{gene} Expression"
        title_para.font.size = Pt(28)
        title_para.font.bold = True
        title_para.alignment = PP_ALIGN.LEFT

        fig_copy = go.Figure(fig)
        fig_copy.update_layout(
            width=1000,
            height=550,
            margin=dict(l=60, r=60, t=60, b=80),
            font=dict(size=14),
        )

        img_bytes = fig_copy.to_image(format="png", scale=2)
        img_stream = io.BytesIO(img_bytes)

        left = Inches(0.5)
        top = Inches(0.9)
        width = Inches(9.5)
        slide.shapes.add_picture(img_stream, left, top, width=width)

        if gene_data is not None and not gene_data.empty:
            stats_box = slide.shapes.add_textbox(
                Inches(10.2), Inches(1.0), Inches(2.8), Inches(5.5)
            )
            stats_frame = stats_box.text_frame
            stats_frame.word_wrap = True

            header_para = stats_frame.paragraphs[0]
            header_para.text = "Statistics"
            header_para.font.size = Pt(14)
            header_para.font.bold = True

            for _, row in gene_data.iterrows():
                cond = row.get("Condition", "N/A")
                fc = row.get("Fold_Change", row.get("Relative_Expression", 0))
                pval = row.get("p_value", float("nan"))
                sig = row.get("significance", "")

                para = stats_frame.add_paragraph()
                para.text = f"{cond[:12]}"
                para.font.size = Pt(10)
                para.font.bold = True

                para2 = stats_frame.add_paragraph()
                fc_str = f"FC: {fc:.2f}" if pd.notna(fc) else "FC: N/A"
                p_str = f"p={pval:.3f}" if pd.notna(pval) else ""
                para2.text = f"  {fc_str} {sig}"
                para2.font.size = Pt(9)

    @staticmethod
    def _add_dual_gene_slide(prs, layout, gene_pairs: list, processed_data: dict):
        from pptx.util import Inches, Pt
        from pptx.enum.text import PP_ALIGN

        slide = prs.slides.add_slide(layout)

        positions = [
            (Inches(0.3), Inches(0.8), Inches(6.2)),
            (Inches(6.8), Inches(0.8), Inches(6.2)),
        ]

        for idx, (gene, fig) in enumerate(gene_pairs):
            if idx >= 2:
                break

            left, top, width = positions[idx]

            title_box = slide.shapes.add_textbox(left, Inches(0.2), width, Inches(0.5))
            title_frame = title_box.text_frame
            title_para = title_frame.paragraphs[0]
            title_para.text = f"{gene}"
            title_para.font.size = Pt(20)
            title_para.font.bold = True

            fig_copy = go.Figure(fig)
            fig_copy.update_layout(
                width=600,
                height=450,
                margin=dict(l=50, r=30, t=40, b=60),
                font=dict(size=11),
            )

            img_bytes = fig_copy.to_image(format="png", scale=2)
            img_stream = io.BytesIO(img_bytes)

            slide.shapes.add_picture(img_stream, left, top, width=width)

    @staticmethod
    def _add_grid_slide(prs, layout, graphs: dict):
        from pptx.util import Inches, Pt
        from pptx.enum.text import PP_ALIGN

        slide = prs.slides.add_slide(layout)

        title_box = slide.shapes.add_textbox(
            Inches(0.3), Inches(0.1), Inches(12.73), Inches(0.5)
        )
        title_frame = title_box.text_frame
        title_para = title_frame.paragraphs[0]
        title_para.text = "Gene Expression Overview"
        title_para.font.size = Pt(24)
        title_para.font.bold = True

        gene_list = list(graphs.items())
        n_genes = len(gene_list)

        if n_genes <= 2:
            cols, rows = 2, 1
        elif n_genes <= 4:
            cols, rows = 2, 2
        elif n_genes <= 6:
            cols, rows = 3, 2
        else:
            cols, rows = 4, 2

        cell_width = 12.5 / cols
        cell_height = 6.5 / rows

        for idx, (gene, fig) in enumerate(gene_list[: cols * rows]):
            col_idx = idx % cols
            row_idx = idx // cols

            left = Inches(0.4 + col_idx * cell_width)
            top = Inches(0.7 + row_idx * cell_height)

            fig_copy = go.Figure(fig)
            fig_copy.update_layout(
                width=350,
                height=280,
                margin=dict(l=40, r=20, t=35, b=40),
                font=dict(size=9),
                title=dict(text=gene, font=dict(size=12)),
            )

            img_bytes = fig_copy.to_image(format="png", scale=2)
            img_stream = io.BytesIO(img_bytes)

            slide.shapes.add_picture(
                img_stream, left, top, width=Inches(cell_width - 0.2)
            )

    @staticmethod
    def _add_summary_slide(prs, layout, processed_data: dict, analysis_params: dict):
        from pptx.util import Inches, Pt
        from pptx.enum.text import PP_ALIGN

        slide = prs.slides.add_slide(layout)

        title_box = slide.shapes.add_textbox(
            Inches(0.3), Inches(0.2), Inches(12.73), Inches(0.6)
        )
        title_frame = title_box.text_frame
        title_para = title_frame.paragraphs[0]
        title_para.text = "Analysis Summary"
        title_para.font.size = Pt(28)
        title_para.font.bold = True

        all_results = (
            pd.concat(processed_data.values(), ignore_index=True)
            if processed_data
            else pd.DataFrame()
        )

        summary_box = slide.shapes.add_textbox(
            Inches(0.5), Inches(1.0), Inches(6), Inches(5.5)
        )
        summary_frame = summary_box.text_frame
        summary_frame.word_wrap = True

        header = summary_frame.paragraphs[0]
        header.text = "Experimental Parameters"
        header.font.size = Pt(16)
        header.font.bold = True

        params_text = [
            f"Efficacy Type: {analysis_params.get('Efficacy_Type', 'N/A')}",
            f"Housekeeping Gene: {analysis_params.get('Housekeeping_Gene', 'N/A')}",
            f"Reference Condition: {analysis_params.get('Reference_Sample', 'N/A')}",
            f"Comparison Condition: {analysis_params.get('Compare_To', 'N/A')}",
            f"Genes Analyzed: {len(processed_data)}",
            f"Analysis Date: {analysis_params.get('Date', 'N/A')}",
        ]

        for text in params_text:
            para = summary_frame.add_paragraph()
            para.text = text
            para.font.size = Pt(12)

        if not all_results.empty:
            results_box = slide.shapes.add_textbox(
                Inches(7), Inches(1.0), Inches(5.8), Inches(5.5)
            )
            results_frame = results_box.text_frame
            results_frame.word_wrap = True

            header2 = results_frame.paragraphs[0]
            header2.text = "Key Findings"
            header2.font.size = Pt(16)
            header2.font.bold = True

            for gene, gene_df in processed_data.items():
                if gene_df.empty:
                    continue

                para = results_frame.add_paragraph()
                para.text = f"\n{gene}:"
                para.font.size = Pt(12)
                para.font.bold = True

                sig_results = gene_df[
                    gene_df.get("significance", pd.Series([""])).str.len() > 0
                ]
                if not sig_results.empty:
                    for _, row in sig_results.head(3).iterrows():
                        cond = row.get("Condition", "N/A")
                        fc = row.get("Fold_Change", row.get("Relative_Expression", 0))
                        sig = row.get("significance", "")
                        para2 = results_frame.add_paragraph()
                        para2.text = f"  {cond}: {fc:.2f}x {sig}"
                        para2.font.size = Pt(10)
                else:
                    para2 = results_frame.add_paragraph()
                    para2.text = "  No significant changes"
                    para2.font.size = Pt(10)


# ==================== EXPORT FUNCTIONS ====================
def export_to_excel(
    raw_data: pd.DataFrame,
    processed_data: Dict[str, pd.DataFrame],
    params: dict,
    mapping: dict,
) -> bytes:
    """Export comprehensive Excel with gene-by-gene sheets"""
    output = io.BytesIO()

    with pd.ExcelWriter(output, engine="xlsxwriter") as writer:
        # Parameters sheet
        pd.DataFrame([params]).to_excel(
            writer, sheet_name="Analysis_Parameters", index=False
        )

        # Sample mapping sheet
        pd.DataFrame([{"Original": k, **v} for k, v in mapping.items()]).to_excel(
            writer, sheet_name="Sample_Mapping", index=False
        )

        # Raw data (include mapped Condition column reflecting sample mapping)
        raw_export = raw_data.copy()
        # mapping provided as argument 'mapping'
        if mapping:
            raw_export["Condition"] = raw_export["Sample"].map(
                lambda x: mapping.get(x, {}).get("condition", x)
            )
        else:
            raw_export["Condition"] = raw_export["Sample"]
        # Keep original Sample name but add mapped Condition
        raw_export = (
            raw_export[["Well", "Sample", "Condition", "Target", "CT", "Source_File"]]
            if "Source_File" in raw_export.columns
            else raw_export
        )
        raw_export.to_excel(writer, sheet_name="Raw_Data", index=False)

        # Gene-by-gene calculations
        for gene, gene_data in processed_data.items():
            sheet_name = f"{gene}_Analysis"[:31]  # Excel sheet name limit
            gene_data.to_excel(writer, sheet_name=sheet_name, index=False)

        # Summary sheet
        if processed_data:
            all_data = pd.concat(processed_data.values(), ignore_index=True)
            summary = (
                all_data.groupby(["Target", "Group"])
                .agg(
                    {"Relative_Expression": ["mean", "std", "count"], "p_value": "min"}
                )
                .round(4)
            )
            summary.to_excel(writer, sheet_name="Summary")

    return output.getvalue()


# ==================== UI ====================
st.title("🧬 qPCR Analysis Suite Pro")
st.markdown("**Gene-by-gene analysis with efficacy-specific workflows**")

# Sidebar
with st.sidebar:
    st.header("💬 Quick Guide")
    st.markdown("""
    1. **📁 Upload** CSV files
    2. **🔍 QC Check** - Review outliers & exclude bad wells
    3. **🗺️ Mapping** - Assign conditions & groups
    4. **🔬 Analysis** - Run ΔΔCt calculations
    5. **📊 Graphs** - Customize visualizations
    6. **📤 Export** - Download publication-ready files
    """)

    st.markdown("---")
    st.markdown("### ⚡ Quick Actions")

    if st.session_state.get("data") is not None:
        st.success(f"✅ {len(st.session_state.data)} wells loaded")

        excluded = len(st.session_state.get("excluded_wells", set()))
        if excluded > 0:
            st.warning(f"⚠️ {excluded} wells excluded")

        if st.session_state.get("processed_data"):
            st.success(f"✅ {len(st.session_state.processed_data)} genes analyzed")
    else:
        st.info("📁 Upload data to begin")

    st.markdown("---")
    with st.expander("⌨️ Navigation Tips"):
        st.markdown("""
        **Tab Navigation**
        - Click tab headers to switch
        - Use Tab/Shift+Tab in forms
        
        **Keyboard Shortcuts**
        - `Ctrl+Enter` - Submit forms
        - `Esc` - Close dialogs
        - `R` - Refresh (browser)
        
        **Quick Tips**
        - Drag column headers to resize tables
        - Double-click graph to reset zoom
        - Hover bars for exact values
        """)

# Main tabs
tab1, tab_qc, tab2, tab3, tab4, tab5, tab6 = st.tabs(
    [
        "📁 Upload",
        "🔍 QC Check",
        "🗺️ Mapping",
        "🔬 Analysis",
        "📊 Graphs",
        "📤 Export",
        "📑 PPT Report",
    ]
)

# ==================== TAB 1: UPLOAD & FILTER ====================
with tab1:
    st.header("Step 1: Upload & Filter Data")

    uploaded_files = st.file_uploader(
        "Upload qPCR CSV files", type=["csv"], accept_multiple_files=True
    )

    if uploaded_files:
        current_file_names = sorted([f.name for f in uploaded_files])
        previous_file_names = st.session_state.get("_uploaded_file_names", [])
        is_new_upload = current_file_names != previous_file_names

        if is_new_upload:
            all_data = []
            for file in uploaded_files:
                parsed = QPCRParser.parse(file)
                if parsed is not None:
                    parsed["Source_File"] = file.name
                    all_data.append(parsed)
                    st.success(f"✅ {file.name}: {len(parsed)} wells")

            if all_data:
                st.session_state.data = pd.concat(all_data, ignore_index=True)
                st.session_state.processed_data = {}
                st.session_state.graphs = {}
                st.session_state.sample_mapping = {}
                st.session_state._uploaded_file_names = current_file_names

                unique_samples = sorted(
                    st.session_state.data["Sample"].unique(), key=natural_sort_key
                )
                st.session_state.sample_order = unique_samples

    if st.session_state.data is not None:
        col1, col2, col3, col4 = st.columns(4)
        col1.metric("Wells", len(st.session_state.data))
        col2.metric("Samples", st.session_state.data["Sample"].nunique())
        col3.metric("Genes", st.session_state.data["Target"].nunique())

        KNOWN_HK_GENES = [
            "ACTIN",
            "B-ACTIN",
            "GAPDH",
            "ACTB",
            "BETA-ACTIN",
            "BETAACTIN",
            "18S",
            "18S RRNA",
            "HPRT",
            "HPRT1",
            "B2M",
            "RPLP0",
            "TBP",
            "PPIA",
            "RPL13A",
            "YWHAZ",
            "SDHA",
            "HMBS",
            "UBC",
            "GUSB",
            "PGK1",
        ]
        all_genes = list(st.session_state.data["Target"].unique())
        hk_genes = [g for g in all_genes if g.upper() in KNOWN_HK_GENES]

        if hk_genes:
            default_idx = 0
            if st.session_state.get("hk_gene") in hk_genes:
                default_idx = hk_genes.index(st.session_state.hk_gene)
            st.session_state.hk_gene = st.selectbox(
                "🔬 Housekeeping Gene (auto-detected)",
                hk_genes,
                index=default_idx,
                key="hk_select",
            )
            col4.metric("HK Gene", st.session_state.hk_gene)
        else:
            st.warning(
                "⚠️ No standard housekeeping gene detected. Please select one manually."
            )
            default_idx = 0
            if st.session_state.get("hk_gene") in all_genes:
                default_idx = all_genes.index(st.session_state.hk_gene)
            st.session_state.hk_gene = st.selectbox(
                "🔬 Select Housekeeping Gene",
                all_genes,
                index=default_idx,
                key="hk_select_manual",
                help="Select the reference/housekeeping gene for normalization",
            )
            col4.metric("HK Gene", st.session_state.hk_gene)

        st.subheader("📊 Data Preview")
        st.dataframe(st.session_state.data.head(50), height=300)

        st.subheader("⚠️ Data Validation")
        warnings_found = False

        data = st.session_state.data
        replicate_counts = data.groupby(["Sample", "Target"]).size()

        single_replicates = replicate_counts[replicate_counts < 2]
        if len(single_replicates) > 0:
            warnings_found = True
            st.warning(
                f"⚠️ **Low Replicates**: {len(single_replicates)} sample-target combinations have only 1 replicate. Statistical analysis requires n≥2."
            )
            with st.expander("View affected samples"):
                st.dataframe(single_replicates.reset_index(name="n"))

        if st.session_state.get("hk_gene"):
            hk = st.session_state.hk_gene
            samples_with_hk = set(data[data["Target"] == hk]["Sample"].unique())
            all_samples = set(data["Sample"].unique())
            missing_hk = all_samples - samples_with_hk
            if missing_hk:
                warnings_found = True
                st.error(
                    f"❌ **Missing Housekeeping Gene**: {len(missing_hk)} samples have no {hk} data: {', '.join(list(missing_hk)[:5])}{'...' if len(missing_hk) > 5 else ''}"
                )

        high_ct_count = len(data[data["CT"] > AnalysisConstants.CT_HIGH_WARNING])
        if high_ct_count > 0:
            warnings_found = True
            pct = high_ct_count / len(data) * 100
            st.warning(
                f"⚠️ **High CT Values**: {high_ct_count} wells ({pct:.1f}%) have CT > {AnalysisConstants.CT_HIGH_WARNING} (low expression)"
            )

        if not warnings_found:
            st.success("✅ Data validation passed. No issues detected.")

# ==================== TAB QC: QUALITY CONTROL ====================
with tab_qc:
    st.header("Step 1.5: Quality Control Check")

    if st.session_state.data is not None and not st.session_state.data.empty:
        data = st.session_state.data
        hk_gene = st.session_state.get("hk_gene")

        if "excluded_wells" not in st.session_state:
            st.session_state.excluded_wells = set()
        if "excluded_wells_history" not in st.session_state:
            st.session_state.excluded_wells_history = []

        # ==================== QC SUMMARY DASHBOARD ====================
        st.markdown(
            """
        <div style='background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); padding: 15px; border-radius: 10px; color: white; margin-bottom: 20px;'>
            <h4 style='margin: 0; color: white;'>🔍 Comprehensive Quality Control</h4>
            <p style='margin: 5px 0 0 0; opacity: 0.9;'>Browse all CT values, review triplicates, and exclude outliers before analysis</p>
        </div>
        """,
            unsafe_allow_html=True,
        )

        # Get comprehensive QC stats
        qc_stats = QualityControl.get_qc_summary_stats(
            data, st.session_state.excluded_wells
        )

        # Summary metrics row
        metric_cols = st.columns(6)
        metric_cols[0].metric("Total Wells", qc_stats.get("total_wells", 0))
        metric_cols[1].metric(
            "Active Wells",
            qc_stats.get("active_wells", 0),
            delta=f"-{qc_stats.get('excluded_wells', 0)}"
            if qc_stats.get("excluded_wells", 0) > 0
            else None,
        )
        metric_cols[2].metric("Triplicates", qc_stats.get("total_triplicates", 0))
        metric_cols[3].metric(
            "Healthy",
            qc_stats.get("healthy_triplicates", 0),
            delta=f"{qc_stats.get('health_score', 0)}%"
            if qc_stats.get("health_score", 0) > 0
            else None,
            delta_color="normal",
        )
        metric_cols[4].metric(
            "⚠️ Warnings",
            qc_stats.get("warning_triplicates", 0),
            delta=None if qc_stats.get("warning_triplicates", 0) == 0 else "review",
            delta_color="off",
        )
        metric_cols[5].metric(
            "❌ Errors",
            qc_stats.get("error_triplicates", 0),
            delta=None if qc_stats.get("error_triplicates", 0) == 0 else "critical",
            delta_color="off",
        )

        # CT Distribution summary
        with st.expander("📈 CT Value Distribution", expanded=False):
            dist_cols = st.columns(4)
            dist_cols[0].metric("Mean CT", f"{qc_stats.get('ct_mean', 0):.1f}")
            dist_cols[1].metric(
                "CT Range",
                f"{qc_stats.get('ct_min', 0):.1f} - {qc_stats.get('ct_max', 0):.1f}",
            )
            dist_cols[2].metric(
                "High CT (>35)",
                qc_stats.get("high_ct_count", 0),
                delta="⚠️" if qc_stats.get("high_ct_count", 0) > 0 else None,
                delta_color="off",
            )
            dist_cols[3].metric(
                "Avg CV%",
                f"{qc_stats.get('avg_cv_pct', 0):.1f}%",
                delta="high" if qc_stats.get("avg_cv_pct", 0) > 5 else None,
                delta_color="off" if qc_stats.get("avg_cv_pct", 0) > 5 else "normal",
            )

        st.markdown("---")

        # ==================== MAIN QC INTERFACE WITH TABS ====================
        qc_tab1, qc_tab2, qc_tab3, qc_tab4 = st.tabs(
            [
                "🔬 Triplicate Browser",
                "🧪 Plate Heatmap",
                "⚠️ Flagged Wells",
                "🔧 Settings",
            ]
        )

        # ==================== TAB 1: TRIPLICATE CT VALUE BROWSER ====================
        with qc_tab1:
            st.subheader("Browse & Edit CT Values by Triplicate")
            st.caption(
                "Scan through all CT values organized by Gene and Sample. Select individual wells to include or exclude from analysis."
            )

            # Filter controls
            filter_col1, filter_col2, filter_col3 = st.columns([1, 1, 1])

            all_genes = sorted(data["Target"].unique())
            all_samples = sorted(data["Sample"].unique(), key=natural_sort_key)

            with filter_col1:
                selected_gene_filter = st.selectbox(
                    "🧬 Filter by Gene", ["All Genes"] + all_genes, key="qc_gene_filter"
                )

            with filter_col2:
                selected_sample_filter = st.selectbox(
                    "🧪 Filter by Sample",
                    ["All Samples"] + all_samples,
                    key="qc_sample_filter",
                )

            with filter_col3:
                status_filter = st.selectbox(
                    "🏷️ Filter by Status",
                    ["All", "OK Only", "Warnings Only", "Errors Only"],
                    key="qc_status_filter",
                )

            # Get triplicate data
            triplicate_data = QualityControl.get_triplicate_data(
                data, st.session_state.excluded_wells
            )

            if not triplicate_data.empty:
                # Apply filters
                filtered_triplicates = triplicate_data.copy()

                if selected_gene_filter != "All Genes":
                    filtered_triplicates = filtered_triplicates[
                        filtered_triplicates["Target"] == selected_gene_filter
                    ]

                if selected_sample_filter != "All Samples":
                    filtered_triplicates = filtered_triplicates[
                        filtered_triplicates["Sample"] == selected_sample_filter
                    ]

                if status_filter == "OK Only":
                    filtered_triplicates = filtered_triplicates[
                        filtered_triplicates["Status"] == "OK"
                    ]
                elif status_filter == "Warnings Only":
                    filtered_triplicates = filtered_triplicates[
                        filtered_triplicates["Severity"] == "warning"
                    ]
                elif status_filter == "Errors Only":
                    filtered_triplicates = filtered_triplicates[
                        filtered_triplicates["Severity"] == "error"
                    ]

                # Show triplicate overview table
                st.markdown(
                    f"**Showing {len(filtered_triplicates)} triplicate groups**"
                )

                # Create a display-friendly dataframe
                display_df = filtered_triplicates[
                    [
                        "Sample",
                        "Target",
                        "n",
                        "Mean_CT",
                        "SD",
                        "CV_pct",
                        "Range",
                        "Status",
                    ]
                ].copy()
                display_df.columns = [
                    "Sample",
                    "Gene",
                    "n",
                    "Mean CT",
                    "SD",
                    "CV%",
                    "Range",
                    "Status",
                ]

                # Style the triplicate overview
                def style_triplicate_row(row):
                    status = row["Status"]
                    if status == "OK":
                        return ["background-color: #d4edda"] * len(row)
                    elif "Has outlier" in status or "High range" in status:
                        return ["background-color: #f8d7da"] * len(row)
                    elif "High CV" in status or "Low n" in status:
                        return ["background-color: #fff3cd"] * len(row)
                    return [""] * len(row)

                styled_overview = display_df.style.apply(
                    style_triplicate_row, axis=1
                ).format(
                    {
                        "Mean CT": "{:.2f}",
                        "SD": "{:.3f}",
                        "CV%": "{:.1f}",
                        "Range": "{:.2f}",
                    }
                )

                st.dataframe(
                    styled_overview,
                    height=250,
                    use_container_width=True,
                    hide_index=True,
                )

                st.markdown("---")

                # Detailed triplicate editor
                st.subheader("📝 Edit Individual Wells Within Triplicates")
                st.caption(
                    "Select a triplicate group to view and edit individual well inclusion/exclusion"
                )

                # Create selection options
                triplicate_options = [
                    f"{row['Sample']} | {row['Target']} (n={row['n']}, CV={row['CV_pct']:.1f}%)"
                    for _, row in filtered_triplicates.iterrows()
                ]

                if triplicate_options:
                    selected_triplicate_idx = st.selectbox(
                        "Select Triplicate Group",
                        range(len(triplicate_options)),
                        format_func=lambda x: triplicate_options[x],
                        key="selected_triplicate",
                    )

                    # Get selected triplicate details
                    selected_row = filtered_triplicates.iloc[selected_triplicate_idx]
                    selected_sample = selected_row["Sample"]
                    selected_target = selected_row["Target"]

                    # Get wells for this triplicate
                    wells_df = QualityControl.get_wells_for_triplicate(
                        data, selected_sample, selected_target
                    )

                    if not wells_df.empty:
                        # Show triplicate summary
                        summary_cols = st.columns([1, 1, 1, 1])
                        summary_cols[0].metric("Gene", selected_target)
                        summary_cols[1].metric(
                            "Sample",
                            selected_sample[:20] + "..."
                            if len(selected_sample) > 20
                            else selected_sample,
                        )
                        summary_cols[2].metric("Status", selected_row["Status"][:30])
                        summary_cols[3].metric(
                            "CV%",
                            f"{selected_row['CV_pct']:.1f}%",
                            delta="⚠️" if selected_row["CV_pct"] > 5 else "✓",
                            delta_color="off"
                            if selected_row["CV_pct"] > 5
                            else "normal",
                        )

                        # Create editable dataframe with Include checkbox
                        editor_df = wells_df.copy()
                        editor_df["Include"] = ~editor_df["Well"].isin(
                            st.session_state.excluded_wells
                        )
                        editor_df["CT"] = editor_df["CT"].round(2)
                        editor_df["Deviation"] = editor_df["Deviation"].round(3)

                        # Reorder columns for display
                        editor_df = editor_df[
                            [
                                "Include",
                                "Well",
                                "CT",
                                "Deviation",
                                "Z_Score",
                                "Well_Status",
                            ]
                        ]

                        # Use data_editor for interactive editing
                        st.markdown(
                            "**Edit well inclusion below** (uncheck 'Include' to exclude a well):"
                        )

                        edited_df = st.data_editor(
                            editor_df,
                            column_config={
                                "Include": st.column_config.CheckboxColumn(
                                    "Include",
                                    help="Uncheck to exclude this well from analysis",
                                    default=True,
                                ),
                                "Well": st.column_config.TextColumn(
                                    "Well",
                                    disabled=True,
                                ),
                                "CT": st.column_config.NumberColumn(
                                    "CT Value",
                                    format="%.2f",
                                    disabled=True,
                                ),
                                "Deviation": st.column_config.NumberColumn(
                                    "Deviation",
                                    help="Deviation from triplicate mean",
                                    format="%.3f",
                                    disabled=True,
                                ),
                                "Z_Score": st.column_config.NumberColumn(
                                    "Z-Score",
                                    help="Standard deviations from mean",
                                    format="%.2f",
                                    disabled=True,
                                ),
                                "Well_Status": st.column_config.TextColumn(
                                    "Status",
                                    disabled=True,
                                ),
                            },
                            hide_index=True,
                            use_container_width=True,
                            key=f"editor_{selected_sample}_{selected_target}",
                        )

                        # Process changes from data_editor
                        if edited_df is not None:
                            for idx, row in edited_df.iterrows():
                                well = row["Well"]
                                include = row["Include"]

                                if (
                                    not include
                                    and well not in st.session_state.excluded_wells
                                ):
                                    # Add to excluded
                                    st.session_state.excluded_wells_history.append(
                                        st.session_state.excluded_wells.copy()
                                    )
                                    st.session_state.excluded_wells.add(well)
                                elif (
                                    include and well in st.session_state.excluded_wells
                                ):
                                    # Remove from excluded
                                    st.session_state.excluded_wells_history.append(
                                        st.session_state.excluded_wells.copy()
                                    )
                                    st.session_state.excluded_wells.discard(well)

                        # Quick action buttons for this triplicate
                        st.markdown("**Quick Actions:**")
                        action_cols = st.columns(4)

                        with action_cols[0]:
                            suggestions = QualityControl.suggest_exclusions(
                                data,
                                selected_sample,
                                selected_target,
                                st.session_state.excluded_wells,
                                strategy="outlier",
                            )
                            if suggestions and st.button(
                                "🎯 Exclude Outliers",
                                key=f"excl_outlier_{selected_sample}_{selected_target}",
                            ):
                                st.session_state.excluded_wells_history.append(
                                    st.session_state.excluded_wells.copy()
                                )
                                for well in suggestions:
                                    st.session_state.excluded_wells.add(well)
                                st.rerun()

                        with action_cols[1]:
                            if st.button(
                                "📉 Keep Best 2",
                                key=f"keep2_{selected_sample}_{selected_target}",
                                help="Keep the 2 CT values closest to median, exclude others",
                            ):
                                suggestions = QualityControl.suggest_exclusions(
                                    data,
                                    selected_sample,
                                    selected_target,
                                    st.session_state.excluded_wells,
                                    strategy="keep_best_2",
                                )
                                if suggestions:
                                    st.session_state.excluded_wells_history.append(
                                        st.session_state.excluded_wells.copy()
                                    )
                                    for well in suggestions:
                                        st.session_state.excluded_wells.add(well)
                                    st.rerun()

                        with action_cols[2]:
                            # Get wells for this triplicate that are currently excluded
                            excluded_in_triplicate = [
                                w
                                for w in wells_df["Well"].tolist()
                                if w in st.session_state.excluded_wells
                            ]
                            if excluded_in_triplicate and st.button(
                                "✅ Include All",
                                key=f"incl_all_{selected_sample}_{selected_target}",
                            ):
                                st.session_state.excluded_wells_history.append(
                                    st.session_state.excluded_wells.copy()
                                )
                                for well in excluded_in_triplicate:
                                    st.session_state.excluded_wells.discard(well)
                                st.rerun()

                        with action_cols[3]:
                            all_wells_in_triplicate = wells_df["Well"].tolist()
                            if st.button(
                                "❌ Exclude All",
                                key=f"excl_all_{selected_sample}_{selected_target}",
                            ):
                                st.session_state.excluded_wells_history.append(
                                    st.session_state.excluded_wells.copy()
                                )
                                for well in all_wells_in_triplicate:
                                    st.session_state.excluded_wells.add(well)
                                st.rerun()
            else:
                st.info("No triplicate data available. Upload data first.")

        # ==================== TAB 2: PLATE HEATMAP ====================
        with qc_tab2:
            st.subheader("96-Well Plate Overview")

            heatmap_col1, heatmap_col2 = st.columns([2, 1])

            with heatmap_col1:
                plate_fig = QualityControl.create_plate_heatmap(
                    data, value_col="CT", excluded_wells=st.session_state.excluded_wells
                )
                st.plotly_chart(plate_fig, use_container_width=True)
                st.caption(
                    "🔴 Red = High CT (low expression) | 🟢 Green = Low CT (high expression) | ❌ = Excluded"
                )

            with heatmap_col2:
                st.markdown("### Replicate Statistics")
                rep_stats = QualityControl.get_replicate_stats(data)
                if not rep_stats.empty:

                    def highlight_status(row):
                        if row["Status"] == "High CV":
                            return ["background-color: #fff3cd"] * len(row)
                        elif row["Status"] == "Low Expression":
                            return ["background-color: #f8d7da"] * len(row)
                        elif row["Status"] == "Check Signal":
                            return ["background-color: #cce5ff"] * len(row)
                        return [""] * len(row)

                    styled_stats = rep_stats.style.apply(highlight_status, axis=1)
                    st.dataframe(styled_stats, height=400, use_container_width=True)

        # ==================== TAB 3: FLAGGED WELLS (LEGACY VIEW) ====================
        with qc_tab3:
            st.subheader("Auto-Flagged Wells")
            st.caption(
                "Wells automatically flagged by QC algorithms (CT thresholds, CV%, Grubbs test)"
            )

            qc_results = QualityControl.detect_outliers(data, hk_gene)
            flagged = qc_results[qc_results["Flagged"]].copy()

            if len(flagged) > 0:
                st.warning(f"Found {len(flagged)} wells with potential issues")

                # Create a more compact display using data_editor
                flagged_display = flagged[
                    ["Well", "Sample", "Target", "CT", "Issues", "Severity"]
                ].copy()
                flagged_display["Exclude"] = flagged_display["Well"].isin(
                    st.session_state.excluded_wells
                )
                flagged_display = flagged_display[
                    ["Exclude", "Well", "Sample", "Target", "CT", "Issues", "Severity"]
                ]

                edited_flagged = st.data_editor(
                    flagged_display,
                    column_config={
                        "Exclude": st.column_config.CheckboxColumn(
                            "Exclude",
                            help="Check to exclude this well from analysis",
                            default=False,
                        ),
                        "CT": st.column_config.NumberColumn(
                            "CT", format="%.2f", disabled=True
                        ),
                        "Well": st.column_config.TextColumn("Well", disabled=True),
                        "Sample": st.column_config.TextColumn("Sample", disabled=True),
                        "Target": st.column_config.TextColumn("Target", disabled=True),
                        "Issues": st.column_config.TextColumn("Issues", disabled=True),
                        "Severity": st.column_config.TextColumn(
                            "Severity", disabled=True
                        ),
                    },
                    hide_index=True,
                    use_container_width=True,
                    key="flagged_wells_editor",
                )

                # Process changes
                if edited_flagged is not None:
                    for idx, row in edited_flagged.iterrows():
                        well = row["Well"]
                        exclude = row["Exclude"]

                        if exclude and well not in st.session_state.excluded_wells:
                            st.session_state.excluded_wells_history.append(
                                st.session_state.excluded_wells.copy()
                            )
                            st.session_state.excluded_wells.add(well)
                        elif not exclude and well in st.session_state.excluded_wells:
                            st.session_state.excluded_wells_history.append(
                                st.session_state.excluded_wells.copy()
                            )
                            st.session_state.excluded_wells.discard(well)

                # Batch action buttons
                batch_cols = st.columns(3)
                with batch_cols[0]:
                    if st.button(
                        "❌ Exclude All Flagged",
                        type="secondary",
                        use_container_width=True,
                    ):
                        st.session_state.excluded_wells_history.append(
                            st.session_state.excluded_wells.copy()
                        )
                        for _, row in flagged.iterrows():
                            st.session_state.excluded_wells.add(row["Well"])
                        st.rerun()

                with batch_cols[1]:
                    if st.button("✅ Clear All Exclusions", use_container_width=True):
                        st.session_state.excluded_wells_history.append(
                            st.session_state.excluded_wells.copy()
                        )
                        st.session_state.excluded_wells = set()
                        st.rerun()

                with batch_cols[2]:
                    can_undo = len(st.session_state.excluded_wells_history) > 0
                    if st.button(
                        "↩️ Undo Last Change",
                        use_container_width=True,
                        disabled=not can_undo,
                    ):
                        if st.session_state.excluded_wells_history:
                            st.session_state.excluded_wells = (
                                st.session_state.excluded_wells_history.pop()
                            )
                            st.rerun()
            else:
                st.success(
                    "✅ No quality issues detected! All wells pass QC thresholds."
                )

        # ==================== TAB 4: QC SETTINGS ====================
        with qc_tab4:
            st.subheader("QC Threshold Settings")
            st.caption("Adjust thresholds for outlier detection algorithms")

            thresh_col1, thresh_col2 = st.columns(2)

            with thresh_col1:
                st.markdown("### CT Value Thresholds")
                new_ct_high = st.number_input(
                    "High CT Threshold",
                    min_value=25.0,
                    max_value=45.0,
                    value=float(QualityControl.CT_HIGH_THRESHOLD),
                    step=0.5,
                    help="Wells with CT above this are flagged as low expression",
                )
                QualityControl.CT_HIGH_THRESHOLD = new_ct_high

                new_ct_low = st.number_input(
                    "Low CT Threshold",
                    min_value=5.0,
                    max_value=20.0,
                    value=float(QualityControl.CT_LOW_THRESHOLD),
                    step=0.5,
                    help="Wells with CT below this are flagged as unusually high expression",
                )
                QualityControl.CT_LOW_THRESHOLD = new_ct_low

            with thresh_col2:
                st.markdown("### Variability Thresholds")
                new_cv = st.number_input(
                    "CV% Threshold",
                    min_value=1.0,
                    max_value=20.0,
                    value=float(QualityControl.CV_THRESHOLD * 100),
                    step=0.5,
                    help="Replicates with CV above this are flagged as high variability",
                )
                QualityControl.CV_THRESHOLD = new_cv / 100

                new_hk_var = st.number_input(
                    "HK Variation Threshold",
                    min_value=0.5,
                    max_value=3.0,
                    value=float(QualityControl.HK_VARIATION_THRESHOLD),
                    step=0.1,
                    help="Housekeeping gene deviation threshold across samples",
                )
                QualityControl.HK_VARIATION_THRESHOLD = new_hk_var

            st.markdown("---")
            st.markdown("### About QC Algorithms")
            st.markdown(
                """
            | Algorithm | Description | Threshold |
            |-----------|-------------|-----------|
            | **CT Range Check** | Flags wells with extreme CT values | CT > {:.0f} (low expression) or CT < {:.0f} (high expression) |
            | **CV% Check** | Flags triplicates with high coefficient of variation | CV > {:.0f}% |
            | **Grubbs Test** | Statistical outlier detection within triplicates (n≥3) | α = {:.2f} |
            | **HK Deviation** | Flags samples where housekeeping gene deviates from mean | Deviation > {:.1f} Ct |
            """.format(
                    QualityControl.CT_HIGH_THRESHOLD,
                    QualityControl.CT_LOW_THRESHOLD,
                    QualityControl.CV_THRESHOLD * 100,
                    QualityControl.GRUBBS_ALPHA,
                    QualityControl.HK_VARIATION_THRESHOLD,
                )
            )

            if st.button(
                "🔄 Re-run QC with New Settings",
                type="primary",
                use_container_width=True,
            ):
                st.rerun()

        # ==================== BOTTOM STATUS BAR ====================
        st.markdown("---")
        excluded_count = len(st.session_state.excluded_wells)

        status_cols = st.columns([2, 1, 1])
        with status_cols[0]:
            if excluded_count > 0:
                st.info(
                    f"ℹ️ **{excluded_count} wells** will be excluded from analysis. Review in Triplicate Browser or proceed to Mapping tab."
                )
            else:
                st.success(
                    "✅ QC check complete. All wells will be included. Proceed to Mapping tab."
                )

        with status_cols[1]:
            can_undo = len(st.session_state.excluded_wells_history) > 0
            if st.button(
                "↩️ Undo Last",
                use_container_width=True,
                disabled=not can_undo,
                key="global_undo",
            ):
                if st.session_state.excluded_wells_history:
                    st.session_state.excluded_wells = (
                        st.session_state.excluded_wells_history.pop()
                    )
                    st.rerun()

        with status_cols[2]:
            if st.button(
                "🔄 Refresh QC", use_container_width=True, key="global_refresh"
            ):
                st.rerun()

    else:
        st.info("⏳ Upload data first in the Upload tab.")

# ==================== TAB 2: SAMPLE MAPPING ====================
with tab2:
    st.header("Step 2: Map Samples to Conditions")

    if st.session_state.data is not None:
        # Efficacy type selection
        detected_genes = set(st.session_state.data["Target"].unique())
        suggested = None
        for eff, cfg in EFFICACY_CONFIG.items():
            if any(g in detected_genes for g in cfg["genes"]):
                suggested = eff
                break

        efficacy = st.selectbox(
            "🎯 Efficacy Test Type",
            list(EFFICACY_CONFIG.keys()),
            index=list(EFFICACY_CONFIG.keys()).index(suggested) if suggested else 0,
        )
        st.session_state.selected_efficacy = efficacy

        config = EFFICACY_CONFIG[efficacy]
        st.info(f"**{config['description']}**")
        st.caption(f"Cell line: {config['cell']} | Genes: {', '.join(config['genes'])}")

        # Show control structure
        with st.expander("📋 Control Structure for this Test"):
            for ctrl_type, ctrl_name in config["controls"].items():
                st.markdown(f"- **{ctrl_type.title()}**: {ctrl_name}")

        # Sample mapping interface with professional layout
        st.markdown("---")
        st.markdown("### 🗺️ Sample Condition Mapping")

        if "sample_order" not in st.session_state or not st.session_state.sample_order:
            st.session_state.sample_order = sorted(
                st.session_state.data["Sample"].unique(), key=natural_sort_key
            )

        # Group type options
        group_types = ["Negative Control", "Positive Control", "Treatment"]
        if "baseline" in config["controls"]:
            group_types.insert(0, "Baseline")

        # Ensure all samples in sample_order have mapping
        for sample in st.session_state.sample_order:
            if sample not in st.session_state.sample_mapping:
                st.session_state.sample_mapping[sample] = {
                    "condition": sample,
                    "group": "Treatment",
                    "concentration": "",
                    "include": True,
                }
            if "include" not in st.session_state.sample_mapping[sample]:
                st.session_state.sample_mapping[sample]["include"] = True

        # Header row with styled background
        st.markdown(
            """
        <div style='background-color: #f8f9fa; padding: 10px; border-radius: 5px; margin-bottom: 10px;'>
            <table style='width: 100%;'>
                <tr>
                    <th style='width: 5%; text-align: center;'>✓</th>
                    <th style='width: 10%;'>Order</th>
                    <th style='width: 15%;'>Original</th>
                    <th style='width: 25%;'>Condition Name</th>
                    <th style='width: 20%;'>Group</th>
                    <th style='width: 10%;'>Move</th>
                </tr>
            </table>
        </div>
        """,
            unsafe_allow_html=True,
        )

        # FIXED: Display ALL samples in sample_order (including excluded ones)
        display_samples = st.session_state.sample_order.copy()

        # Sample rows with improved spacing
        for i, sample in enumerate(display_samples):
            # Container for each row
            with st.container():
                col0, col_order, col1, col2, col3, col_move = st.columns(
                    [0.5, 0.8, 1.5, 2.5, 2, 1]
                )

                # Include checkbox
                with col0:
                    include = st.checkbox(
                        "Include sample",
                        value=st.session_state.sample_mapping[sample].get(
                            "include", True
                        ),
                        key=f"include_{sample}_{i}",
                        label_visibility="collapsed",
                    )
                    st.session_state.sample_mapping[sample]["include"] = include

                # Order number
                with col_order:
                    st.markdown(
                        f"<div style='text-align: center; padding-top: 10px;'><b>{i + 1}</b></div>",
                        unsafe_allow_html=True,
                    )

                # Original sample name (non-editable)
                with col1:
                    st.text_input(
                        "Original",
                        sample,
                        key=f"orig_{sample}_{i}",
                        disabled=True,
                        label_visibility="collapsed",
                    )

                # Condition name (editable)
                with col2:
                    cond = st.text_input(
                        "Condition",
                        st.session_state.sample_mapping[sample]["condition"],
                        key=f"cond_{sample}_{i}",
                        label_visibility="collapsed",
                        placeholder="Enter condition name...",
                    )
                    st.session_state.sample_mapping[sample]["condition"] = cond

                # Group selector
                with col3:
                    grp_idx = 0
                    try:
                        grp_idx = group_types.index(
                            st.session_state.sample_mapping[sample]["group"]
                        )
                    except (ValueError, KeyError):
                        pass

                    grp = st.selectbox(
                        "Group",
                        group_types,
                        index=grp_idx,
                        key=f"grp_{sample}_{i}",
                        label_visibility="collapsed",
                    )
                    st.session_state.sample_mapping[sample]["group"] = grp

                # Move controls - FIXED: Use immutable operations to prevent race conditions
                with col_move:
                    btn_col1, btn_col2 = st.columns(2)
                    with btn_col1:
                        if i > 0:
                            if st.button(
                                "⬆",
                                key=f"up_{sample}_{i}",
                                help="Move up",
                                width="stretch",
                            ):
                                # Create new list with swapped items (immutable operation)
                                new_order = st.session_state.sample_order.copy()
                                new_order[i], new_order[i - 1] = (
                                    new_order[i - 1],
                                    new_order[i],
                                )
                                st.session_state.sample_order = new_order
                                st.rerun()
                    with btn_col2:
                        if i < len(display_samples) - 1:
                            if st.button(
                                "⬇",
                                key=f"down_{sample}_{i}",
                                help="Move down",
                                width="stretch",
                            ):
                                # Create new list with swapped items (immutable operation)
                                new_order = st.session_state.sample_order.copy()
                                new_order[i], new_order[i + 1] = (
                                    new_order[i + 1],
                                    new_order[i],
                                )
                                st.session_state.sample_order = new_order
                                st.rerun()

                # Divider line
                st.markdown(
                    "<hr style='margin: 5px 0; opacity: 0.3;'>", unsafe_allow_html=True
                )

        # Update excluded_samples from include flags
        st.session_state.excluded_samples = set(
            [
                s
                for s, v in st.session_state.sample_mapping.items()
                if not v.get("include", True)
            ]
        )

        # Summary with styled cards
        st.markdown("---")
        st.subheader("📊 Mapping Summary")

        col_card1, col_card2, col_card3, col_card4 = st.columns(4)

        total_samples = len(st.session_state.sample_order)
        included = sum(
            1
            for s in st.session_state.sample_order
            if st.session_state.sample_mapping[s].get("include", True)
        )
        excluded = total_samples - included

        with col_card1:
            st.metric("Total Samples", total_samples)
        with col_card2:
            st.metric(
                "Included",
                included,
                delta=None if included == total_samples else f"-{excluded}",
            )
        with col_card3:
            st.metric(
                "Excluded", excluded, delta=None if excluded == 0 else f"+{excluded}"
            )
        with col_card4:
            groups = set(
                v["group"]
                for v in st.session_state.sample_mapping.values()
                if v.get("include", True)
            )
            st.metric("Groups", len(groups))

        # Detailed table view
        with st.expander("📋 View Detailed Mapping Table"):
            mapping_df = pd.DataFrame(
                [
                    {
                        "Order": idx + 1,
                        "Include": "✅"
                        if st.session_state.sample_mapping[s].get("include", True)
                        else "❌",
                        "Original": s,
                        "Condition": st.session_state.sample_mapping[s]["condition"],
                        "Group": st.session_state.sample_mapping[s]["group"],
                    }
                    for idx, s in enumerate(st.session_state.sample_order)
                ]
            )
            st.dataframe(mapping_df, width="stretch", hide_index=True)

        # Run analysis
        st.markdown("---")
        st.subheader("🔬 Run Full Analysis (ΔΔCt + Statistics)")

        # Build condition list from mapping
        condition_list = []
        sample_to_condition = {}

        for sample in st.session_state.get("sample_order", []):
            if sample in st.session_state.sample_mapping:
                mapping_info = st.session_state.sample_mapping[sample]
                if mapping_info.get("include", True):
                    condition = mapping_info.get("condition", sample)
                    condition_list.append(condition)
                    sample_to_condition[condition] = sample

        if condition_list:
            # Enhanced layout with clear separation
            st.markdown("#### 📊 Analysis Configuration")

            col_info1, col_info2 = st.columns(2)
            with col_info1:
                st.info(
                    "**ΔΔCt Reference:** Used to calculate fold changes. All samples will be relative to this (Fold Change = 1.0)"
                )
            with col_info2:
                st.info(
                    "**P-value References:** Used for statistical comparison (t-test). Choose one or two conditions for comparison."
                )

            col_r1, col_r2, col_r3 = st.columns(3)
            with col_r1:
                ref_condition = st.selectbox(
                    "🎯 ΔΔCt Reference Condition",
                    condition_list,
                    index=0,
                    key="ref_choice_ddct",
                    help="Baseline for relative expression calculation",
                )
                ref_sample_key = sample_to_condition[ref_condition]
                st.caption(f"→ Sample: **{ref_sample_key}**")

            with col_r2:
                cmp_condition = st.selectbox(
                    "📈 P-value Reference 1 (*)",
                    condition_list,
                    index=0,
                    key="cmp_choice_pval",
                    help="Primary control group for statistical testing (asterisk symbols)",
                )
                cmp_sample_key = sample_to_condition[cmp_condition]
                st.caption(f"→ Sample: **{cmp_sample_key}**")

            with col_r3:
                # Add option for second p-value comparison
                use_second_comparison = st.checkbox(
                    "Enable 2nd comparison (#)",
                    value=False,
                    key="use_second_pval",
                    help="Add a second statistical comparison with hashtag symbols",
                )

                if use_second_comparison:
                    # Filter out the first comparison from options
                    condition_list_2 = [c for c in condition_list if c != cmp_condition]
                    if condition_list_2:
                        cmp_condition_2 = st.selectbox(
                            "📊 P-value Reference 2 (#)",
                            condition_list_2,
                            index=0,
                            key="cmp_choice_pval_2",
                            help="Secondary control group for statistical testing (hashtag symbols)",
                        )
                        cmp_sample_key_2 = sample_to_condition[cmp_condition_2]
                        st.caption(f"→ Sample: **{cmp_sample_key_2}**")
                    else:
                        st.warning("Need at least 3 conditions for dual comparison")
                        use_second_comparison = False
                        cmp_sample_key_2 = None
                else:
                    cmp_sample_key_2 = None

            # Visual summary
            st.markdown("---")
            col_sum1, col_sum2, col_sum3 = st.columns([1, 2, 1])
            with col_sum2:
                summary_html = f"""
                <div style='background-color: #f0f2f6; padding: 15px; border-radius: 10px; text-align: center;'>
                    <h4>Analysis Summary</h4>
                    <p><b>Fold Changes:</b> Relative to <code>{ref_condition}</code></p>
                    <p><b>P-values (*):</b> Compared to <code>{cmp_condition}</code></p>
                """
                if use_second_comparison and cmp_sample_key_2:
                    summary_html += f"<p><b>P-values (#):</b> Compared to <code>{cmp_condition_2}</code></p>"
                summary_html += "</div>"
                st.markdown(summary_html, unsafe_allow_html=True)

            st.markdown("---")

            # Run button
            if st.button("▶️ Run Full Analysis Now", type="primary", width="stretch"):
                ok = AnalysisEngine.run_full_analysis(
                    ref_sample_key,
                    cmp_sample_key,
                    cmp_sample_key_2 if use_second_comparison else None,
                )
                if ok:
                    success_msg = f"✅ Analysis complete!\n\n- Fold changes relative to: **{ref_condition}**\n- P-values (*) vs: **{cmp_condition}**"
                    if use_second_comparison and cmp_sample_key_2:
                        success_msg += f"\n- P-values (#) vs: **{cmp_condition_2}**"
                    st.success(success_msg)
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
        all_results = pd.concat(
            st.session_state.processed_data.values(), ignore_index=True
        )

        col1, col2, col3 = st.columns(3)
        col1.metric("Genes Analyzed", len(st.session_state.processed_data))
        col2.metric("Conditions", all_results["Condition"].nunique())
        sig_count = (all_results["p_value"] < 0.05).sum()
        col3.metric("Significant (p<0.05)", f"{sig_count}/{len(all_results)}")

        # Show results per gene
        st.subheader("🧬 Gene-by-Gene Results")

        for gene, gene_df in st.session_state.processed_data.items():
            with st.expander(f"📍 {gene}", expanded=False):
                # Show expected direction if available
                efficacy_config = EFFICACY_CONFIG.get(
                    st.session_state.selected_efficacy, {}
                )
                if "expected_direction" in efficacy_config:
                    direction = efficacy_config["expected_direction"].get(gene)
                    if direction:
                        st.caption(
                            f"Expected: {'↑ Increase' if direction == 'up' else '↓ Decrease'}"
                        )

                # Display columns
                display_cols = [
                    "Condition",
                    "Group",
                    "Fold_Change",
                    "p_value",
                    "significance",
                    "n_replicates",
                    "Target_Ct_Mean",
                    "HK_Ct_Mean",
                    "Delta_Ct",
                    "SEM",
                ]

                # Filter to existing columns
                display_df = gene_df[[c for c in display_cols if c in gene_df.columns]]

                # Style the dataframe
                styled = display_df.style.background_gradient(
                    subset=["Fold_Change"], cmap="RdYlGn", vmin=0, vmax=3
                ).format(
                    {
                        "Fold_Change": "{:.3f}",
                        "p_value": "{:.4f}",
                        "Target_Ct_Mean": "{:.2f}",
                        "HK_Ct_Mean": "{:.2f}",
                        "Delta_Ct": "{:.2f}",
                        "SEM": "{:.3f}",
                    },
                    na_rep="—",
                )

                st.dataframe(styled, width="stretch")

        st.success("✅ Results ready! Go to Graphs tab to visualize.")

    else:
        st.info(
            "⏳ No analysis results yet. Go to 'Sample Mapping' tab and click 'Run Full Analysis Now'"
        )

# ==================== TAB 4: GRAPHS ====================
with tab4:
    st.header("Step 4: Individual Gene Graphs")

    st.markdown(
        """
    <style>
    [data-testid="stPlotlyChart"] {
        box-shadow: 0 4px 12px rgba(0,0,0,0.15);
        border-radius: 8px;
        background: white;
        padding: 10px;
    }
    
    .stExpander {
        background-color: #FAFAFA;
        border-left: 3px solid #E0E0E0;
        margin-bottom: 5px;
    }
    
    .gene-selector-btn {
        transition: all 0.2s ease;
    }
    
    .graph-toolbar {
        background: linear-gradient(135deg, #f5f7fa 0%, #e8ecf1 100%);
        padding: 12px 16px;
        border-radius: 8px;
        margin-bottom: 16px;
    }
    
    .stat-highlight {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        padding: 4px 8px;
        border-radius: 4px;
        font-weight: 600;
    }
    </style>
    """,
        unsafe_allow_html=True,
    )

    if st.session_state.processed_data:
        if "graph_settings" not in st.session_state:
            st.session_state.graph_settings = {
                "title_size": 20,
                "font_size": 14,
                "sig_font_size": 18,
                "figure_width": 1000,
                "figure_height": 600,
                "color_scheme": "plotly_white",
                "show_error": True,
                "show_significance": True,
                "show_grid": True,
                "xlabel": "Condition",
                "ylabel": "Relative mRNA Expression Level",
                "bar_colors": {},
                "orientation": "v",
                "error_multiplier": 1.96,
                "bar_opacity": 0.95,
                "bar_gap": 0.15,
                "marker_line_width": 1,
                "show_legend": False,
                "y_log_scale": False,
                "y_min": None,
                "y_max": None,
            }

        if "graphs" not in st.session_state:
            st.session_state.graphs = {}

        efficacy_config = EFFICACY_CONFIG.get(st.session_state.selected_efficacy, {})
        gene_list = list(st.session_state.processed_data.keys())

        if "selected_gene_idx" not in st.session_state:
            st.session_state.selected_gene_idx = 0

        st.markdown(
            """
        <style>
        .gene-pill { display: inline-block; padding: 8px 16px; margin: 2px; border-radius: 20px;
                     font-weight: 600; cursor: pointer; transition: all 0.2s; }
        .gene-pill-active { background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; }
        .gene-pill-inactive { background: #f0f2f6; color: #333; }
        .compact-control { background: #fafafa; padding: 8px; border-radius: 8px; margin: 4px 0; }
        .color-grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(80px, 1fr)); gap: 4px; }
        </style>
        """,
            unsafe_allow_html=True,
        )

        st.markdown("### 🧬 Select Gene")
        gene_cols = st.columns(min(len(gene_list), 6))
        for idx, gene in enumerate(gene_list):
            with gene_cols[idx % len(gene_cols)]:
                if st.button(
                    f"{'✓ ' if idx == st.session_state.selected_gene_idx else ''}{gene}",
                    key=f"gene_btn_{gene}",
                    width="stretch",
                    type="primary"
                    if idx == st.session_state.selected_gene_idx
                    else "secondary",
                ):
                    st.session_state.selected_gene_idx = idx
                    st.rerun()

        current_gene = gene_list[st.session_state.selected_gene_idx]
        gene_data = st.session_state.processed_data[current_gene]

        st.markdown("---")

        toolbar_cols = st.columns([1, 1, 1, 1, 2])

        show_sig_key = f"{current_gene}_show_sig"
        show_err_key = f"{current_gene}_show_err"
        bar_gap_key = f"{current_gene}_bar_gap"

        if show_sig_key not in st.session_state.graph_settings:
            st.session_state.graph_settings[show_sig_key] = True
        if show_err_key not in st.session_state.graph_settings:
            st.session_state.graph_settings[show_err_key] = True
        if bar_gap_key not in st.session_state.graph_settings:
            st.session_state.graph_settings[bar_gap_key] = 0.25

        with toolbar_cols[0]:
            sig_on = st.toggle(
                "✨ Significance",
                st.session_state.graph_settings[show_sig_key],
                key=f"tgl_sig_{current_gene}",
            )
            st.session_state.graph_settings[show_sig_key] = sig_on

        with toolbar_cols[1]:
            err_on = st.toggle(
                "📏 Error Bars",
                st.session_state.graph_settings[show_err_key],
                key=f"tgl_err_{current_gene}",
            )
            st.session_state.graph_settings[show_err_key] = err_on

        with toolbar_cols[2]:
            gap_val = st.select_slider(
                "Gap",
                options=[0.1, 0.15, 0.2, 0.25, 0.3, 0.4],
                value=st.session_state.graph_settings[bar_gap_key],
                key=f"gap_sl_{current_gene}",
            )
            st.session_state.graph_settings[bar_gap_key] = gap_val

        with toolbar_cols[3]:
            if st.button(
                "↺ Reset All", key=f"reset_all_{current_gene}", width="stretch"
            ):
                if f"{current_gene}_bar_settings" in st.session_state:
                    del st.session_state[f"{current_gene}_bar_settings"]
                st.session_state.graph_settings[show_sig_key] = True
                st.session_state.graph_settings[show_err_key] = True
                st.session_state.graph_settings[bar_gap_key] = 0.25
                st.rerun()

        with toolbar_cols[4]:
            edit_mode = st.checkbox(
                "🎨 Edit Bar Colors", key=f"edit_mode_{current_gene}"
            )

        if f"{current_gene}_bar_settings" not in st.session_state:
            st.session_state[f"{current_gene}_bar_settings"] = {}
        if "bar_colors_per_sample" not in st.session_state.graph_settings:
            st.session_state.graph_settings["bar_colors_per_sample"] = {}

        for idx, (_, row) in enumerate(gene_data.iterrows()):
            condition = row["Condition"]
            group = row.get("Group", "Treatment")
            default_color = DEFAULT_GROUP_COLORS.get(group, "#D3D3D3")
            bar_key = f"{current_gene}_{condition}"

            if bar_key not in st.session_state[f"{current_gene}_bar_settings"]:
                st.session_state[f"{current_gene}_bar_settings"][bar_key] = {
                    "color": default_color,
                    "show_sig": True,
                    "show_err": True,
                }

        if edit_mode:
            with st.expander("🎨 Bar Color & Visibility Editor", expanded=True):
                st.markdown(
                    """
                <style>
                .color-editor-card {
                    background: #f8f9fa;
                    border: 1px solid #e9ecef;
                    border-radius: 8px;
                    padding: 12px;
                    margin-bottom: 8px;
                }
                .color-editor-title {
                    font-size: 13px;
                    font-weight: 600;
                    color: #333;
                    margin-bottom: 4px;
                    white-space: nowrap;
                    overflow: hidden;
                    text-overflow: ellipsis;
                }
                .color-editor-subtitle {
                    font-size: 11px;
                    color: #666;
                    margin-bottom: 8px;
                }
                </style>
                """,
                    unsafe_allow_html=True,
                )

                n_bars = len(gene_data)
                n_cols = min(n_bars, 3)

                st.markdown("#### Individual Bar Settings")
                color_cols = st.columns(n_cols)

                for idx, (_, row) in enumerate(gene_data.iterrows()):
                    condition = row["Condition"]
                    group = row.get("Group", "Treatment")
                    bar_key = f"{current_gene}_{condition}"

                    with color_cols[idx % n_cols]:
                        with st.container():
                            display_name = (
                                condition
                                if len(condition) <= 18
                                else condition[:15] + "..."
                            )
                            st.markdown(
                                f"""
                            <div class='color-editor-card'>
                                <div class='color-editor-title' title='{condition}'>{display_name}</div>
                                <div class='color-editor-subtitle'>{group}</div>
                            </div>
                            """,
                                unsafe_allow_html=True,
                            )

                            new_color = st.color_picker(
                                "Color",
                                st.session_state[f"{current_gene}_bar_settings"][
                                    bar_key
                                ]["color"],
                                key=f"cp_{current_gene}_{idx}",
                                help=f"Bar color for {condition}",
                            )
                            st.session_state[f"{current_gene}_bar_settings"][bar_key][
                                "color"
                            ] = new_color
                            st.session_state.graph_settings["bar_colors_per_sample"][
                                bar_key
                            ] = new_color

                            cb_col1, cb_col2 = st.columns(2)
                            with cb_col1:
                                sig_bar = st.checkbox(
                                    "Show *",
                                    st.session_state[f"{current_gene}_bar_settings"][
                                        bar_key
                                    ]["show_sig"],
                                    key=f"sb_{current_gene}_{idx}",
                                    help="Show significance marker",
                                )
                                st.session_state[f"{current_gene}_bar_settings"][
                                    bar_key
                                ]["show_sig"] = sig_bar
                            with cb_col2:
                                err_bar = st.checkbox(
                                    "Show ±",
                                    st.session_state[f"{current_gene}_bar_settings"][
                                        bar_key
                                    ]["show_err"],
                                    key=f"eb_{current_gene}_{idx}",
                                    help="Show error bar",
                                )
                                st.session_state[f"{current_gene}_bar_settings"][
                                    bar_key
                                ]["show_err"] = err_bar

                            st.markdown(
                                "<div style='margin-bottom: 16px;'></div>",
                                unsafe_allow_html=True,
                            )

                st.markdown("---")
                st.markdown("#### Quick Color Presets")
                st.caption("Apply a color scheme to all bars")

                preset_cols = st.columns(4)
                with preset_cols[0]:
                    if st.button(
                        "🔵 Blues",
                        key=f"preset_blue_{current_gene}",
                        use_container_width=True,
                    ):
                        blues = ["#e3f2fd", "#90caf9", "#42a5f5", "#1976d2", "#0d47a1"]
                        for idx, (_, row) in enumerate(gene_data.iterrows()):
                            bar_key = f"{current_gene}_{row['Condition']}"
                            st.session_state[f"{current_gene}_bar_settings"][bar_key][
                                "color"
                            ] = blues[idx % len(blues)]
                            st.session_state.graph_settings["bar_colors_per_sample"][
                                bar_key
                            ] = blues[idx % len(blues)]
                        st.rerun()
                with preset_cols[1]:
                    if st.button(
                        "🟢 Greens",
                        key=f"preset_green_{current_gene}",
                        use_container_width=True,
                    ):
                        greens = ["#e8f5e9", "#a5d6a7", "#66bb6a", "#388e3c", "#1b5e20"]
                        for idx, (_, row) in enumerate(gene_data.iterrows()):
                            bar_key = f"{current_gene}_{row['Condition']}"
                            st.session_state[f"{current_gene}_bar_settings"][bar_key][
                                "color"
                            ] = greens[idx % len(greens)]
                            st.session_state.graph_settings["bar_colors_per_sample"][
                                bar_key
                            ] = greens[idx % len(greens)]
                        st.rerun()
                with preset_cols[2]:
                    if st.button(
                        "🔴 Warm",
                        key=f"preset_warm_{current_gene}",
                        use_container_width=True,
                    ):
                        warms = ["#fff3e0", "#ffcc80", "#ff9800", "#f57c00", "#e65100"]
                        for idx, (_, row) in enumerate(gene_data.iterrows()):
                            bar_key = f"{current_gene}_{row['Condition']}"
                            st.session_state[f"{current_gene}_bar_settings"][bar_key][
                                "color"
                            ] = warms[idx % len(warms)]
                            st.session_state.graph_settings["bar_colors_per_sample"][
                                bar_key
                            ] = warms[idx % len(warms)]
                        st.rerun()
                with preset_cols[3]:
                    if st.button(
                        "⬜ Grayscale",
                        key=f"preset_gray_{current_gene}",
                        use_container_width=True,
                    ):
                        grays = ["#ffffff", "#e0e0e0", "#bdbdbd", "#9e9e9e", "#616161"]
                        for idx, (_, row) in enumerate(gene_data.iterrows()):
                            bar_key = f"{current_gene}_{row['Condition']}"
                            st.session_state[f"{current_gene}_bar_settings"][bar_key][
                                "color"
                            ] = grays[idx % len(grays)]
                            st.session_state.graph_settings["bar_colors_per_sample"][
                                bar_key
                            ] = grays[idx % len(grays)]
                        st.rerun()

        current_settings = st.session_state.graph_settings.copy()
        current_settings["show_significance"] = st.session_state.graph_settings.get(
            show_sig_key, True
        )
        current_settings["show_error"] = st.session_state.graph_settings.get(
            show_err_key, True
        )
        current_settings["bar_gap"] = st.session_state.graph_settings.get(
            bar_gap_key, 0.25
        )

        fig = GraphGenerator.create_gene_graph(
            gene_data,
            current_gene,
            current_settings,
            efficacy_config,
            sample_order=st.session_state.get("sample_order"),
            per_sample_overrides=None,
        )

        st.plotly_chart(fig, width="stretch", key=f"main_fig_{current_gene}")
        st.session_state.graphs[current_gene] = fig

        with st.expander("📊 All Gene Graphs (Quick View)", expanded=False):
            all_gene_cols = st.columns(min(len(gene_list), 2))
            for idx, gene in enumerate(gene_list):
                if gene == current_gene:
                    continue

                gd = st.session_state.processed_data[gene]

                if f"{gene}_bar_settings" not in st.session_state:
                    st.session_state[f"{gene}_bar_settings"] = {}
                    for _, row in gd.iterrows():
                        condition = row["Condition"]
                        group = row.get("Group", "Treatment")
                        default_color = DEFAULT_GROUP_COLORS.get(group, "#D3D3D3")
                        bar_key = f"{gene}_{condition}"
                        st.session_state[f"{gene}_bar_settings"][bar_key] = {
                            "color": default_color,
                            "show_sig": True,
                            "show_err": True,
                        }

                gs = st.session_state.graph_settings.copy()
                gs["show_significance"] = gs.get(f"{gene}_show_sig", True)
                gs["show_error"] = gs.get(f"{gene}_show_err", True)
                gs["bar_gap"] = gs.get(f"{gene}_bar_gap", 0.25)
                gs["figure_height"] = 350

                f = GraphGenerator.create_gene_graph(
                    gd,
                    gene,
                    gs,
                    efficacy_config,
                    sample_order=st.session_state.get("sample_order"),
                )

                with all_gene_cols[idx % len(all_gene_cols)]:
                    st.markdown(f"**{gene}**")
                    st.plotly_chart(f, width="stretch", key=f"mini_{gene}")
                    st.session_state.graphs[gene] = f
    else:
        st.info(
            "⏳ No analysis results yet. Go to 'Sample Mapping' tab and click 'Run Full Analysis Now'"
        )

# ==================== TAB 5: EXPORT ====================
with tab5:
    st.header("Step 5: Export All Results")

    if st.session_state.processed_data:
        st.subheader("📦 Download Options")

        analysis_params = {
            "Date": datetime.now().strftime("%Y-%m-%d %H:%M"),
            "Efficacy_Type": st.session_state.selected_efficacy,
            "Housekeeping_Gene": st.session_state.hk_gene,
            "Reference_Sample": st.session_state.get("analysis_ref_condition", "N/A"),
            "Compare_To": st.session_state.get("analysis_cmp_condition", "N/A"),
            "Excluded_Wells": len(st.session_state.excluded_wells),
            "Excluded_Samples": len(st.session_state.excluded_samples),
            "Genes_Analyzed": len(st.session_state.processed_data),
        }

        col1, col2 = st.columns(2)

        with col1:
            st.markdown("### 📊 Complete Excel Report")
            st.caption(
                "Includes: Parameters, Mapping, Raw Data, Gene-by-Gene Calculations, Summary"
            )

            excel_data = export_to_excel(
                st.session_state.data,
                st.session_state.processed_data,
                analysis_params,
                st.session_state.sample_mapping,
            )

            st.download_button(
                label="📥 Download Excel Report",
                data=excel_data,
                file_name=f"qPCR_{st.session_state.selected_efficacy}_{datetime.now().strftime('%Y%m%d_%H%M')}.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                type="primary",
            )

        with col2:
            st.markdown("### 📈 All Graphs (HTML)")
            st.caption("Interactive graphs for all genes in one file")

            if st.session_state.graphs:
                # Create combined HTML
                html_parts = [
                    "<html><head><title>qPCR Analysis Graphs</title></head><body>"
                ]
                html_parts.append(
                    f"<h1>{st.session_state.selected_efficacy} Analysis</h1>"
                )
                html_parts.append(
                    f"<p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>"
                )

                for gene, fig in st.session_state.graphs.items():
                    html_parts.append(f"<h2>{gene}</h2>")
                    html_parts.append(
                        fig.to_html(include_plotlyjs="cdn", div_id=f"graph_{gene}")
                    )
                    html_parts.append("<hr>")

                html_parts.append("</body></html>")
                combined_html = "\n".join(html_parts)

                st.download_button(
                    label="📥 Download All Graphs (HTML)",
                    data=combined_html,
                    file_name=f"qPCR_graphs_{st.session_state.selected_efficacy}_{datetime.now().strftime('%Y%m%d_%H%M')}.html",
                    mime="text/html",
                    type="primary",
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
                    key=f"csv_{gene}",
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
                    key=f"html_{gene}",
                )

        with col3:
            # Configuration export
            st.markdown("**Reproducibility Files**")

            # Analysis config
            config_data = {
                "analysis_params": analysis_params,
                "sample_mapping": st.session_state.sample_mapping,
                "graph_settings": st.session_state.graph_settings,
                "excluded_wells": list(st.session_state.excluded_wells),
                "excluded_samples": list(st.session_state.excluded_samples),
            }

            st.download_button(
                label="📥 Analysis Config",
                data=json.dumps(config_data, indent=2),
                file_name=f"config_{datetime.now().strftime('%Y%m%d_%H%M')}.json",
                mime="application/json",
            )

            # Graph preset
            st.download_button(
                label="📥 Graph Preset",
                data=json.dumps(st.session_state.graph_settings, indent=2),
                file_name=f"graph_preset_{datetime.now().strftime('%Y%m%d')}.json",
                mime="application/json",
            )

        st.markdown("---")
        st.subheader("📸 Publication-Ready Images")
        st.caption("High-resolution images suitable for journals and reports")

        pub_col1, pub_col2, pub_col3 = st.columns(3)

        with pub_col1:
            img_format = st.selectbox(
                "Image Format",
                ["PNG (300 DPI)", "SVG (Vector)", "PDF (Vector)"],
                key="pub_img_format",
            )

        with pub_col2:
            img_width = st.number_input(
                "Width (px)", min_value=400, max_value=3000, value=1200, step=100
            )

        with pub_col3:
            img_height = st.number_input(
                "Height (px)", min_value=300, max_value=2000, value=800, step=100
            )

        if st.session_state.graphs:
            st.markdown("**Download Individual High-Res Images:**")

            img_cols = st.columns(min(len(st.session_state.graphs), 4))

            for idx, (gene, fig) in enumerate(st.session_state.graphs.items()):
                col_idx = idx % len(img_cols)
                with img_cols[col_idx]:
                    try:
                        fig_copy = go.Figure(fig)
                        fig_copy.update_layout(
                            width=img_width,
                            height=img_height,
                            font=dict(size=14),
                            title=dict(font=dict(size=18)),
                        )

                        if "PNG" in img_format:
                            img_bytes = fig_copy.to_image(
                                format="png",
                                scale=3,
                                width=img_width,
                                height=img_height,
                            )
                            st.download_button(
                                label=f"📥 {gene}.png",
                                data=img_bytes,
                                file_name=f"{gene}_{datetime.now().strftime('%Y%m%d')}_300dpi.png",
                                mime="image/png",
                                key=f"png_{gene}",
                            )
                        elif "SVG" in img_format:
                            svg_bytes = fig_copy.to_image(
                                format="svg", width=img_width, height=img_height
                            )
                            st.download_button(
                                label=f"📥 {gene}.svg",
                                data=svg_bytes,
                                file_name=f"{gene}_{datetime.now().strftime('%Y%m%d')}.svg",
                                mime="image/svg+xml",
                                key=f"svg_{gene}",
                            )
                        else:
                            pdf_bytes = fig_copy.to_image(
                                format="pdf", width=img_width, height=img_height
                            )
                            st.download_button(
                                label=f"📥 {gene}.pdf",
                                data=pdf_bytes,
                                file_name=f"{gene}_{datetime.now().strftime('%Y%m%d')}.pdf",
                                mime="application/pdf",
                                key=f"pdf_{gene}",
                            )
                    except Exception as e:
                        st.warning(
                            f"Image export requires kaleido: pip install kaleido"
                        )
                        break

            st.markdown("---")
            st.subheader("📦 Batch Export (ZIP)")

            batch_col1, batch_col2 = st.columns(2)

            with batch_col1:
                if st.button(
                    "📥 Download All Figures (ZIP)", type="primary", width="stretch"
                ):
                    try:
                        zip_buffer = io.BytesIO()
                        total_graphs = len(st.session_state.graphs)
                        progress_bar = st.progress(0, text="Preparing figures...")

                        with zipfile.ZipFile(
                            zip_buffer, "w", zipfile.ZIP_DEFLATED
                        ) as zf:
                            for i, (gene, fig) in enumerate(
                                st.session_state.graphs.items()
                            ):
                                progress_bar.progress(
                                    (i + 1) / total_graphs,
                                    text=f"Exporting {gene}... ({i + 1}/{total_graphs})",
                                )
                                fig_copy = go.Figure(fig)
                                fig_copy.update_layout(
                                    width=img_width, height=img_height
                                )

                                png_bytes = fig_copy.to_image(
                                    format="png",
                                    scale=3,
                                    width=img_width,
                                    height=img_height,
                                )
                                zf.writestr(f"{gene}_300dpi.png", png_bytes)

                                svg_bytes = fig_copy.to_image(
                                    format="svg", width=img_width, height=img_height
                                )
                                zf.writestr(f"{gene}.svg", svg_bytes)

                                html_str = fig_copy.to_html(include_plotlyjs="cdn")
                                zf.writestr(f"{gene}.html", html_str)

                        progress_bar.empty()
                        st.download_button(
                            label="📥 Download ZIP Now",
                            data=zip_buffer.getvalue(),
                            file_name=f"qPCR_figures_{st.session_state.selected_efficacy}_{datetime.now().strftime('%Y%m%d')}.zip",
                            mime="application/zip",
                            key="batch_zip",
                        )
                        st.success(
                            f"✅ Created ZIP with {len(st.session_state.graphs)} figures (PNG + SVG + HTML)"
                        )
                    except Exception as e:
                        st.error(f"Batch export requires kaleido: pip install kaleido")

            with batch_col2:
                if st.button("📥 Download Complete Report (ZIP)", width="stretch"):
                    try:
                        zip_buffer = io.BytesIO()
                        total_steps = (
                            2
                            + len(st.session_state.graphs)
                            + len(st.session_state.processed_data)
                        )
                        current_step = 0
                        progress_bar = st.progress(0, text="Creating report...")

                        with zipfile.ZipFile(
                            zip_buffer, "w", zipfile.ZIP_DEFLATED
                        ) as zf:
                            progress_bar.progress(
                                1 / total_steps, text="Generating Excel report..."
                            )
                            current_step += 1
                            zf.writestr(
                                "analysis_report.xlsx",
                                export_to_excel(
                                    st.session_state.data,
                                    st.session_state.processed_data,
                                    analysis_params,
                                    st.session_state.sample_mapping,
                                ),
                            )

                            progress_bar.progress(
                                2 / total_steps, text="Saving config..."
                            )
                            current_step += 1
                            zf.writestr(
                                "analysis_config.json",
                                json.dumps(
                                    {
                                        "analysis_params": analysis_params,
                                        "sample_mapping": st.session_state.sample_mapping,
                                        "excluded_wells": list(
                                            st.session_state.excluded_wells
                                        ),
                                        "excluded_samples": list(
                                            st.session_state.excluded_samples
                                        ),
                                    },
                                    indent=2,
                                ),
                            )

                            for i, (gene, fig) in enumerate(
                                st.session_state.graphs.items()
                            ):
                                current_step += 1
                                progress_bar.progress(
                                    current_step / total_steps,
                                    text=f"Exporting figure: {gene}...",
                                )
                                fig_copy = go.Figure(fig)
                                fig_copy.update_layout(
                                    width=img_width, height=img_height
                                )

                                png_bytes = fig_copy.to_image(
                                    format="png",
                                    scale=3,
                                    width=img_width,
                                    height=img_height,
                                )
                                zf.writestr(f"figures/{gene}_300dpi.png", png_bytes)

                                svg_bytes = fig_copy.to_image(
                                    format="svg", width=img_width, height=img_height
                                )
                                zf.writestr(f"figures/{gene}.svg", svg_bytes)

                            for (
                                gene,
                                gene_df,
                            ) in st.session_state.processed_data.items():
                                current_step += 1
                                progress_bar.progress(
                                    current_step / total_steps,
                                    text=f"Exporting data: {gene}...",
                                )
                                csv_str = gene_df.to_csv(index=False)
                                zf.writestr(f"data/{gene}_results.csv", csv_str)

                        progress_bar.empty()
                        st.download_button(
                            label="📥 Download Complete Report ZIP",
                            data=zip_buffer.getvalue(),
                            file_name=f"qPCR_complete_{st.session_state.selected_efficacy}_{datetime.now().strftime('%Y%m%d')}.zip",
                            mime="application/zip",
                            key="complete_zip",
                        )
                        st.success("✅ Complete report package created!")
                    except Exception as e:
                        st.error(f"Report generation failed: {e}")

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

# ==================== TAB 6: PPT REPORT ====================
with tab6:
    st.header("Step 6: PowerPoint Presentation Export")

    pptx_available = False
    try:
        from pptx import Presentation as _TestPptx

        pptx_available = True
    except ImportError:
        pass

    if not pptx_available:
        st.error(
            "⚠️ PowerPoint export is temporarily unavailable. The python-pptx package failed to load."
        )
        with st.expander("🔧 Debug Info"):
            import sys

            st.code(f"Python version: {sys.version}")
            try:
                import pkg_resources

                installed = [f"{d.key}=={d.version}" for d in pkg_resources.working_set]
                pptx_pkgs = [
                    p for p in installed if "pptx" in p.lower() or "lxml" in p.lower()
                ]
                st.write("**Relevant packages:**")
                if pptx_pkgs:
                    for p in pptx_pkgs:
                        st.code(p)
                else:
                    st.warning("python-pptx and lxml NOT found in installed packages")
            except Exception as e:
                st.error(f"Could not check packages: {e}")
    else:
        st.success("✅ PowerPoint export ready")

    st.markdown(
        """
    <div style='background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); padding: 20px; border-radius: 12px; color: white; margin-bottom: 20px;'>
        <h3 style='margin: 0; color: white;'>📑 Publication-Ready Presentations</h3>
        <p style='margin: 8px 0 0 0; opacity: 0.9;'>Generate professional PowerPoint slides with your gene expression graphs</p>
    </div>
    """,
        unsafe_allow_html=True,
    )

    if st.session_state.graphs and st.session_state.processed_data and pptx_available:
        st.subheader("🎨 Presentation Settings")

        col_layout, col_options = st.columns(2)

        with col_layout:
            st.markdown("#### Slide Layout")
            layout_option = st.radio(
                "Choose how to arrange graphs:",
                ["one_per_slide", "two_per_slide", "grid"],
                format_func=lambda x: {
                    "one_per_slide": "📊 One graph per slide (best for presentations)",
                    "two_per_slide": "📊📊 Two graphs per slide (compact)",
                    "grid": "📋 All graphs on one slide (overview)",
                }[x],
                key="ppt_layout",
            )

            st.markdown("#### Slide Content")
            include_title = st.checkbox(
                "Include title slide", value=True, key="ppt_title"
            )
            include_summary = st.checkbox(
                "Include summary slide", value=True, key="ppt_summary"
            )

        with col_options:
            st.markdown("#### Preview")

            n_genes = len(st.session_state.graphs)

            if layout_option == "one_per_slide":
                n_slides = (
                    n_genes
                    + (1 if include_title else 0)
                    + (1 if include_summary else 0)
                )
                st.info(
                    f"📄 **{n_slides} slides** will be generated:\n"
                    f"- {'1 title slide' if include_title else ''}\n"
                    f"- {n_genes} gene slides\n"
                    f"- {'1 summary slide' if include_summary else ''}"
                )
            elif layout_option == "two_per_slide":
                gene_slides = (n_genes + 1) // 2
                n_slides = (
                    gene_slides
                    + (1 if include_title else 0)
                    + (1 if include_summary else 0)
                )
                st.info(
                    f"📄 **{n_slides} slides** will be generated:\n"
                    f"- {'1 title slide' if include_title else ''}\n"
                    f"- {gene_slides} gene slides (2 per slide)\n"
                    f"- {'1 summary slide' if include_summary else ''}"
                )
            else:
                n_slides = (
                    1 + (1 if include_title else 0) + (1 if include_summary else 0)
                )
                st.info(
                    f"📄 **{n_slides} slides** will be generated:\n"
                    f"- {'1 title slide' if include_title else ''}\n"
                    f"- 1 overview slide (all genes)\n"
                    f"- {'1 summary slide' if include_summary else ''}"
                )

            st.markdown("#### Genes to Include")
            gene_list = list(st.session_state.graphs.keys())
            st.write(", ".join([f"**{g}**" for g in gene_list]))

        st.markdown("---")

        st.subheader("📥 Generate & Download")

        gen_col1, gen_col2, gen_col3 = st.columns([2, 1, 1])

        with gen_col1:
            if st.button(
                "🚀 Generate PowerPoint Presentation",
                type="primary",
                use_container_width=True,
            ):
                try:
                    with st.spinner("Generating presentation..."):
                        analysis_params = {
                            "Date": datetime.now().strftime("%Y-%m-%d %H:%M"),
                            "Efficacy_Type": st.session_state.get(
                                "selected_efficacy", "qPCR Analysis"
                            ),
                            "Housekeeping_Gene": st.session_state.get("hk_gene", "N/A"),
                            "Reference_Sample": st.session_state.get(
                                "analysis_ref_condition", "N/A"
                            ),
                            "Compare_To": st.session_state.get(
                                "analysis_cmp_condition", "N/A"
                            ),
                            "Genes_Analyzed": len(st.session_state.processed_data),
                        }

                        ppt_bytes = ReportGenerator.create_presentation(
                            graphs=st.session_state.graphs,
                            processed_data=st.session_state.processed_data,
                            analysis_params=analysis_params,
                            layout=layout_option,
                            include_title_slide=include_title,
                            include_summary=include_summary,
                        )

                        st.session_state["ppt_bytes"] = ppt_bytes
                        st.success("✅ Presentation generated successfully!")

                except ImportError as e:
                    st.error(f"❌ Import error: {str(e)}")
                    st.info(
                        "💡 If on Streamlit Cloud, try rebooting the app to reinstall dependencies."
                    )
                except Exception as e:
                    st.error(f"❌ Error generating presentation: {str(e)}")

        if "ppt_bytes" in st.session_state and st.session_state["ppt_bytes"]:
            with gen_col2:
                efficacy = st.session_state.get("selected_efficacy", "qPCR")
                filename = (
                    f"qPCR_{efficacy}_{datetime.now().strftime('%Y%m%d_%H%M')}.pptx"
                )

                st.download_button(
                    label="📥 Download PPTX",
                    data=st.session_state["ppt_bytes"],
                    file_name=filename,
                    mime="application/vnd.openxmlformats-officedocument.presentationml.presentation",
                    type="primary",
                    use_container_width=True,
                )

        st.markdown("---")

        with st.expander("💡 Tips for Great Presentations", expanded=False):
            st.markdown("""
            **Layout Recommendations:**
            - **One per slide**: Best for formal presentations, allows detailed discussion of each gene
            - **Two per slide**: Good for comparing related genes or showing more data in less time
            - **Grid view**: Great for overview slides or quick reference handouts
            
            **After Download:**
            1. Open in PowerPoint or Google Slides
            2. Adjust fonts and colors to match your organization's template
            3. Add additional context or notes as needed
            4. Consider adding a methods slide manually for complete presentations
            
            **For Publications:**
            - Individual high-resolution images can be downloaded from the Export tab
            - SVG format is recommended for vector graphics in publications
            """)

    elif pptx_available:
        st.info(
            "⏳ Generate graphs first in the Graphs tab to create a PowerPoint presentation."
        )

        st.markdown("""
        ### How to use:
        1. **Upload** your qPCR data in the Upload tab
        2. **Map** your samples to conditions in the Mapping tab
        3. **Run Analysis** to calculate ΔΔCt values
        4. **Generate Graphs** in the Graphs tab
        5. **Return here** to create your PowerPoint presentation
        """)

# ==================== FOOTER ====================
st.markdown("---")

footer_html = """
<div style='text-align: center; color: #666;'>
    <p>🧬 qPCR Analysis Suite Pro v3.1 | Gene-by-gene analysis with efficacy-specific workflows</p>
    <p>QC Check • Outlier Detection • Publication-Ready Export • PowerPoint Reports</p>
</div>
"""

st.markdown(footer_html, unsafe_allow_html=True)
