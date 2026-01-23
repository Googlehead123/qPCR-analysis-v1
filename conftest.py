"""
Pytest configuration and fixtures for qPCR Analysis Platform tests.

This module provides shared fixtures and mocks for testing the qPCR analysis
application without requiring Streamlit runtime.
"""

import sys
from unittest.mock import MagicMock

import pytest
import pandas as pd
import numpy as np


# ==================== STREAMLIT MOCK ====================
# Mock streamlit before importing the main module
class MockSessionState(dict):
    """Mock Streamlit session_state that behaves like a dict with attribute access."""

    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError:
            raise AttributeError(f"'MockSessionState' object has no attribute '{key}'")

    def __setattr__(self, key, value):
        self[key] = value

    def __delattr__(self, key):
        try:
            del self[key]
        except KeyError:
            raise AttributeError(f"'MockSessionState' object has no attribute '{key}'")


class MockContextManager:
    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass


def _create_mock_streamlit():
    mock_st = MagicMock()
    mock_st.session_state = MockSessionState()
    mock_st.error = MagicMock()
    mock_st.warning = MagicMock()
    mock_st.success = MagicMock()
    mock_st.info = MagicMock()
    mock_st.spinner = MagicMock(return_value=MockContextManager())
    mock_st.tabs = MagicMock(return_value=[MockContextManager() for _ in range(7)])
    mock_st.sidebar = MockContextManager()
    mock_st.columns = MagicMock(return_value=[MagicMock() for _ in range(3)])
    mock_st.expander = MagicMock(return_value=MockContextManager())
    mock_st.form = MagicMock(return_value=MockContextManager())
    mock_st.set_page_config = MagicMock()
    mock_st.title = MagicMock()
    mock_st.markdown = MagicMock()
    mock_st.header = MagicMock()
    mock_st.subheader = MagicMock()
    mock_st.write = MagicMock()
    mock_st.dataframe = MagicMock()
    mock_st.plotly_chart = MagicMock()
    mock_st.file_uploader = MagicMock(return_value=None)
    mock_st.selectbox = MagicMock(return_value=None)
    mock_st.multiselect = MagicMock(return_value=[])
    mock_st.text_input = MagicMock(return_value="")
    mock_st.number_input = MagicMock(return_value=0)
    mock_st.checkbox = MagicMock(return_value=False)
    mock_st.button = MagicMock(return_value=False)
    mock_st.radio = MagicMock(return_value=None)
    mock_st.slider = MagicMock(return_value=0)
    mock_st.color_picker = MagicMock(return_value="#FFFFFF")
    mock_st.download_button = MagicMock(return_value=False)
    mock_st.metric = MagicMock()
    mock_st.rerun = MagicMock()
    mock_st.cache_data = lambda f: f
    mock_st.toggle = MagicMock(return_value=True)
    mock_st.data_editor = MagicMock(return_value=None)
    mock_st.caption = MagicMock()
    mock_st.select_slider = MagicMock(return_value=0)
    mock_st.column_config = MagicMock()
    mock_st.container = MagicMock(return_value=MockContextManager())
    return mock_st


sys.modules["streamlit"] = _create_mock_streamlit()


@pytest.fixture(autouse=True)
def mock_streamlit():
    """Auto-use fixture to mock Streamlit for all tests."""
    mock_st = _create_mock_streamlit()
    sys.modules["streamlit"] = mock_st

    main_module_name = "streamlit qpcr analysis v1"
    if main_module_name in sys.modules:
        del sys.modules[main_module_name]

    yield mock_st


# ==================== SAMPLE DATA FIXTURES ====================
@pytest.fixture
def sample_qpcr_raw_data():
    """Generate realistic qPCR raw data for testing.

    Returns a DataFrame mimicking parsed qPCR CSV output with:
    - 3 samples (Non-treated, Treatment1, Treatment2)
    - 3 replicates each
    - 2 targets (GAPDH housekeeping, COL1A1 gene of interest)
    """
    data = []

    samples = ["Non-treated", "Treatment1", "Treatment2"]
    targets = ["GAPDH", "COL1A1"]

    # Base Ct values (realistic qPCR values)
    base_ct = {
        ("Non-treated", "GAPDH"): 18.5,
        ("Non-treated", "COL1A1"): 25.0,
        ("Treatment1", "GAPDH"): 18.3,
        ("Treatment1", "COL1A1"): 23.5,  # ~3x upregulation vs Non-treated
        ("Treatment2", "GAPDH"): 18.6,
        ("Treatment2", "COL1A1"): 26.5,  # ~0.35x downregulation
    }

    well_counter = 1
    for sample in samples:
        for target in targets:
            for rep in range(3):
                ct_value = base_ct[(sample, target)] + np.random.normal(0, 0.2)
                well = f"A{well_counter}"
                well_counter += 1

                data.append(
                    {
                        "Well": well,
                        "Sample": sample,
                        "Target": target,
                        "CT": round(ct_value, 2),
                    }
                )

    return pd.DataFrame(data)


@pytest.fixture
def sample_mapping():
    """Sample mapping configuration for tests."""
    return {
        "Non-treated": {
            "condition": "Non-treated",
            "group": "Negative Control",
            "include": True,
        },
        "Treatment1": {
            "condition": "Treatment1",
            "group": "Treatment",
            "include": True,
        },
        "Treatment2": {
            "condition": "Treatment2",
            "group": "Treatment",
            "include": True,
        },
    }


@pytest.fixture
def format1_csv_content():
    """CSV content in Format 1 (standard QuantStudio export)."""
    return """Block Type,,,,,,,
Experiment File Name,test.eds,,,,,,
Experiment Run End Time,2024-01-15 10:30:00,,,,,,

Well Position,Well,Sample Name,Target Name,Task,Reporter,Quencher,CT
A1,A1,Non-treated,GAPDH,UNKNOWN,FAM,NFQ-MGB,18.52
A2,A2,Non-treated,GAPDH,UNKNOWN,FAM,NFQ-MGB,18.48
A3,A3,Non-treated,GAPDH,UNKNOWN,FAM,NFQ-MGB,18.55
A4,A4,Non-treated,COL1A1,UNKNOWN,FAM,NFQ-MGB,25.12
A5,A5,Non-treated,COL1A1,UNKNOWN,FAM,NFQ-MGB,24.89
A6,A6,Non-treated,COL1A1,UNKNOWN,FAM,NFQ-MGB,25.01
B1,B1,Treatment1,GAPDH,UNKNOWN,FAM,NFQ-MGB,18.33
B2,B2,Treatment1,GAPDH,UNKNOWN,FAM,NFQ-MGB,18.28
B3,B3,Treatment1,GAPDH,UNKNOWN,FAM,NFQ-MGB,18.35
B4,B4,Treatment1,COL1A1,UNKNOWN,FAM,NFQ-MGB,23.45
B5,B5,Treatment1,COL1A1,UNKNOWN,FAM,NFQ-MGB,23.52
B6,B6,Treatment1,COL1A1,UNKNOWN,FAM,NFQ-MGB,23.48
"""


@pytest.fixture
def format2_csv_content():
    """CSV content in Format 2 (alternative export with Cyrillic Ct symbol)."""
    return """Well,Well Position,Omit,Sample Name,Target Name,Task,Reporter,Quencher,Ct,Ct Mean,Ct SD,Quantity,Quantity Mean,Quantity SD,Y-Intercept,R^2,Slope,Efficiency,Automatic Ct Threshold,Ct Threshold,Automatic Baseline,Baseline Start,Baseline End,Amp Status,Comments,HIGHSD,NOAMP,OUTLIERRG,EXPFAIL,Tm1,Tm2,Tm3
A1,A1,FALSE,Non-treated,GAPDH,UNKNOWN,FAM,NFQ-MGB,18.52,,,,,,,,,,TRUE,0.2,TRUE,3,15,Amp,,FALSE,FALSE,FALSE,FALSE,,,
A2,A2,FALSE,Non-treated,GAPDH,UNKNOWN,FAM,NFQ-MGB,18.48,,,,,,,,,,TRUE,0.2,TRUE,3,15,Amp,,FALSE,FALSE,FALSE,FALSE,,,
A3,A3,FALSE,Non-treated,COL1A1,UNKNOWN,FAM,NFQ-MGB,25.12,,,,,,,,,,TRUE,0.2,TRUE,3,15,Amp,,FALSE,FALSE,FALSE,FALSE,,,
"""


@pytest.fixture
def malformed_csv_content():
    """Malformed CSV with insufficient columns."""
    return """Col1,Col2,Col3
A1,Sample1,Value1
A2,Sample2,Value2
"""


@pytest.fixture
def processed_gene_data():
    """Pre-processed gene data as would come from AnalysisEngine.calculate_ddct()."""
    return pd.DataFrame(
        {
            "Target": ["COL1A1", "COL1A1", "COL1A1"],
            "Condition": ["Non-treated", "Treatment1", "Treatment2"],
            "Original_Sample": ["Non-treated", "Treatment1", "Treatment2"],
            "Group": ["Negative Control", "Treatment", "Treatment"],
            "n_replicates": [3, 3, 3],
            "Target_Ct_Mean": [25.01, 23.48, 26.50],
            "Target_Ct_SD": [0.12, 0.04, 0.15],
            "HK_Ct_Mean": [18.52, 18.32, 18.57],
            "Delta_Ct": [6.49, 5.16, 7.93],
            "Delta_Delta_Ct": [0.0, -1.33, 1.44],
            "Relative_Expression": [1.0, 2.52, 0.37],
            "SEM": [0.07, 0.02, 0.09],
            "Fold_Change": [1.0, 2.52, 0.37],
        }
    )


@pytest.fixture
def graph_settings():
    """Default graph settings for visualization tests."""
    return {
        "show_error": True,
        "show_significance": True,
        "error_multiplier": 1.96,
        "bar_colors": {"COL1A1": "#4CAF50"},
        "bar_colors_per_sample": {},
        "marker_line_width": 1,
        "bar_opacity": 0.95,
        "bar_gap": 0.15,
        "title_size": 20,
        "font_size": 14,
        "figure_height": 600,
        "figure_width": 1000,
        "color_scheme": "plotly_white",
        "plot_bgcolor": "#FFFFFF",
        "show_legend": False,
        "y_log_scale": False,
        "y_min": None,
        "y_max": None,
    }
