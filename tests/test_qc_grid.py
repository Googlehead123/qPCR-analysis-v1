"""
RED phase unit tests for QC grid state management functions.

These tests are designed to FAIL initially (TDD approach) because the
functions being tested don't exist yet. They define the expected behavior
for grid cell selection state management.

Grid state design:
- Session state key: 'qc_grid_selected_cell'
- Value: {'gene': str, 'sample': str} or None
- Functions to implement:
  - get_grid_cell_key(gene, sample) -> str
  - get_selected_cell(session_state) -> tuple or None
  - set_selected_cell(session_state, gene, sample) -> None
  - clear_selected_cell(session_state) -> None
  - is_cell_selected(session_state, gene, sample) -> bool
"""

import pytest
from importlib import import_module


class TestQCGridStateManagement:
    """Test suite for QC grid cell selection state management functions."""

    def test_get_grid_cell_key_creates_unique_key(self, mock_streamlit):
        """
        get_grid_cell_key should create a unique, consistent string key
        from gene and sample names.

        This key is used internally to identify grid cells and ensure
        consistent state management across UI updates.
        """
        spec = import_module("streamlit qpcr analysis v1")
        get_grid_cell_key = spec.get_grid_cell_key

        # Same inputs should produce same key
        key1 = get_grid_cell_key("GENE1", "Sample_A")
        key2 = get_grid_cell_key("GENE1", "Sample_A")
        assert key1 == key2, "Same inputs should produce identical keys"

        # Different inputs should produce different keys
        key3 = get_grid_cell_key("GENE2", "Sample_A")
        key4 = get_grid_cell_key("GENE1", "Sample_B")
        assert key1 != key3, "Different genes should produce different keys"
        assert key1 != key4, "Different samples should produce different keys"

        # Key should be a string
        assert isinstance(key1, str), "Key must be a string"
        assert len(key1) > 0, "Key must not be empty"

    def test_set_and_get_selected_cell(self, mock_streamlit):
        """
        set_selected_cell should store the selected cell in session state,
        and get_selected_cell should retrieve it correctly.

        This validates the basic read/write cycle for grid selection state.
        """
        spec = import_module("streamlit qpcr analysis v1")
        set_selected_cell = spec.set_selected_cell
        get_selected_cell = spec.get_selected_cell

        session_state = mock_streamlit.session_state

        # Initially, no cell should be selected
        result = get_selected_cell(session_state)
        assert result is None, "No cell should be selected initially"

        # Set a cell as selected
        set_selected_cell(session_state, "GENE1", "Sample_A")

        # Retrieve the selected cell
        selected = get_selected_cell(session_state)
        assert selected is not None, "Selected cell should not be None"
        assert selected == ("GENE1", "Sample_A"), (
            "Selected cell should return (gene, sample) tuple"
        )

    def test_set_selected_cell_overwrites_previous(self, mock_streamlit):
        """
        set_selected_cell should overwrite the previously selected cell.

        Only one cell can be selected at a time. Setting a new cell
        should replace the old selection.
        """
        spec = import_module("streamlit qpcr analysis v1")
        set_selected_cell = spec.set_selected_cell
        get_selected_cell = spec.get_selected_cell

        session_state = mock_streamlit.session_state

        # Set first cell
        set_selected_cell(session_state, "GENE1", "Sample_A")
        assert get_selected_cell(session_state) == ("GENE1", "Sample_A")

        # Set different cell - should overwrite
        set_selected_cell(session_state, "GENE2", "Sample_B")
        selected = get_selected_cell(session_state)
        assert selected == ("GENE2", "Sample_B"), (
            "New selection should overwrite previous selection"
        )

    def test_clear_selected_cell(self, mock_streamlit):
        """
        clear_selected_cell should remove the current selection,
        returning the grid to an unselected state.

        This is used when user clicks away or when filters change.
        """
        spec = import_module("streamlit qpcr analysis v1")
        set_selected_cell = spec.set_selected_cell
        clear_selected_cell = spec.clear_selected_cell
        get_selected_cell = spec.get_selected_cell

        session_state = mock_streamlit.session_state

        # Set a cell
        set_selected_cell(session_state, "GENE1", "Sample_A")
        assert get_selected_cell(session_state) is not None

        # Clear the selection
        clear_selected_cell(session_state)

        # Verify it's cleared
        result = get_selected_cell(session_state)
        assert result is None, "Selected cell should be None after clearing"

    def test_is_cell_selected_true_when_matching(self, mock_streamlit):
        """
        is_cell_selected should return True when the given gene and sample
        match the currently selected cell.

        This is used to highlight the selected cell in the grid UI.
        """
        spec = import_module("streamlit qpcr analysis v1")
        set_selected_cell = spec.set_selected_cell
        is_cell_selected = spec.is_cell_selected

        session_state = mock_streamlit.session_state

        # Set a cell as selected
        set_selected_cell(session_state, "GENE1", "Sample_A")

        # Check if the same cell is selected
        assert is_cell_selected(session_state, "GENE1", "Sample_A") is True, (
            "Should return True for the selected cell"
        )

    def test_is_cell_selected_false_when_not_matching(self, mock_streamlit):
        """
        is_cell_selected should return False when the given gene and sample
        do NOT match the currently selected cell.

        This ensures only the correct cell is highlighted in the grid.
        """
        spec = import_module("streamlit qpcr analysis v1")
        set_selected_cell = spec.set_selected_cell
        is_cell_selected = spec.is_cell_selected

        session_state = mock_streamlit.session_state

        # Set a cell as selected
        set_selected_cell(session_state, "GENE1", "Sample_A")

        # Check different cells
        assert is_cell_selected(session_state, "GENE2", "Sample_A") is False, (
            "Should return False for different gene"
        )
        assert is_cell_selected(session_state, "GENE1", "Sample_B") is False, (
            "Should return False for different sample"
        )
        assert is_cell_selected(session_state, "GENE2", "Sample_B") is False, (
            "Should return False for completely different cell"
        )

    def test_is_cell_selected_false_when_nothing_selected(self, mock_streamlit):
        """
        is_cell_selected should return False when no cell is currently selected.

        This handles the initial state where the grid has no selection.
        """
        spec = import_module("streamlit qpcr analysis v1")
        is_cell_selected = spec.is_cell_selected

        session_state = mock_streamlit.session_state

        # No cell has been selected yet
        assert is_cell_selected(session_state, "GENE1", "Sample_A") is False, (
            "Should return False when nothing is selected"
        )

    def test_cell_selection_independence(self, mock_streamlit):
        """
        Cell selection should be independent - selecting one cell should NOT
        affect the selection state of other cells.

        This validates that the state management doesn't have side effects
        that could cause unexpected behavior when multiple cells are involved.
        """
        spec = import_module("streamlit qpcr analysis v1")
        set_selected_cell = spec.set_selected_cell
        is_cell_selected = spec.is_cell_selected

        session_state = mock_streamlit.session_state

        # Select cell A
        set_selected_cell(session_state, "GENE1", "Sample_A")
        assert is_cell_selected(session_state, "GENE1", "Sample_A") is True

        # Verify cell B is NOT selected
        assert is_cell_selected(session_state, "GENE1", "Sample_B") is False, (
            "Selecting cell A should not affect cell B selection state"
        )

        # Select cell B
        set_selected_cell(session_state, "GENE1", "Sample_B")
        assert is_cell_selected(session_state, "GENE1", "Sample_B") is True

        # Verify cell A is now NOT selected (overwritten)
        assert is_cell_selected(session_state, "GENE1", "Sample_A") is False, (
            "Selecting cell B should overwrite cell A selection"
        )

    def test_state_persistence_across_calls(self, mock_streamlit):
        """
        State should persist across multiple function calls.

        This validates that the session_state is properly maintained
        and not reset between operations.
        """
        spec = import_module("streamlit qpcr analysis v1")
        set_selected_cell = spec.set_selected_cell
        get_selected_cell = spec.get_selected_cell
        is_cell_selected = spec.is_cell_selected

        session_state = mock_streamlit.session_state

        # Set a cell
        set_selected_cell(session_state, "GENE1", "Sample_A")

        # Make multiple calls - state should persist
        for _ in range(3):
            assert get_selected_cell(session_state) == ("GENE1", "Sample_A")
            assert is_cell_selected(session_state, "GENE1", "Sample_A") is True

        # State should still be there
        assert get_selected_cell(session_state) == ("GENE1", "Sample_A")


class TestQCGridUIHelpers:
    """Test suite for QC grid UI helper functions.

    These tests validate data transformation functions that convert
    triplicate-level qPCR data into grid display structures.

    Helper functions to implement:
    - build_grid_matrix(triplicate_data) -> dict
    - get_cell_status_color(status) -> str
    - get_cell_display_text(cell_data) -> str
    """

    def test_build_grid_matrix_creates_nested_dict_structure(
        self, mock_streamlit, sample_qpcr_raw_data
    ):
        """
        build_grid_matrix should transform triplicate DataFrame into
        a nested dictionary structure: {gene: {sample: cell_data}}.

        This structure enables efficient grid rendering where each cell
        contains aggregated statistics for a gene-sample combination.
        """
        spec = import_module("streamlit qpcr analysis v1")
        build_grid_matrix = spec.build_grid_matrix

        # Get triplicate data (requires QualityControl class)
        qc = spec.QualityControl()
        triplicate_data = qc.get_triplicate_data(sample_qpcr_raw_data)

        # Build grid matrix
        matrix = build_grid_matrix(triplicate_data)

        # Should return a dictionary
        assert isinstance(matrix, dict), "Grid matrix must be a dictionary"

        # Should have genes as top-level keys
        assert len(matrix) > 0, "Grid matrix should not be empty"

        # Each gene should map to a dict of samples
        for gene, samples_dict in matrix.items():
            assert isinstance(gene, str), "Gene key must be a string"
            assert isinstance(samples_dict, dict), (
                "Each gene should map to a dictionary of samples"
            )

            # Each sample should have cell data
            for sample, cell_data in samples_dict.items():
                assert isinstance(sample, str), "Sample key must be a string"
                assert isinstance(cell_data, dict), "Cell data must be a dictionary"

    def test_build_grid_matrix_cell_data_contains_required_fields(
        self, mock_streamlit, sample_qpcr_raw_data
    ):
        """
        Each cell in the grid matrix should contain required fields:
        mean_ct, cv, status, n.

        These fields are necessary for rendering cell content and
        determining cell styling.
        """
        spec = import_module("streamlit qpcr analysis v1")
        build_grid_matrix = spec.build_grid_matrix

        qc = spec.QualityControl()
        triplicate_data = qc.get_triplicate_data(sample_qpcr_raw_data)

        matrix = build_grid_matrix(triplicate_data)

        # Check that all cells have required fields
        required_fields = {"mean_ct", "cv", "status", "n"}
        for gene, samples_dict in matrix.items():
            for sample, cell_data in samples_dict.items():
                for field in required_fields:
                    assert field in cell_data, (
                        f"Cell data for {gene}/{sample} missing field: {field}"
                    )

    def test_get_cell_status_color_maps_ok_to_green(self, mock_streamlit):
        """
        get_cell_status_color should map "OK" status to green color.

        Green (#d4edda) indicates the cell has no quality issues.
        """
        spec = import_module("streamlit qpcr analysis v1")
        get_cell_status_color = spec.get_cell_status_color

        color = get_cell_status_color("OK")

        assert isinstance(color, str), "Color must be a string"
        assert color == "#d4edda", "OK status should map to green color (#d4edda)"

    def test_get_cell_status_color_maps_warnings_to_yellow(self, mock_streamlit):
        """
        get_cell_status_color should map warning statuses to yellow color.

        Yellow (#fff3cd) indicates quality concerns that need attention
        but are not critical (e.g., High CV, Low n).
        """
        spec = import_module("streamlit qpcr analysis v1")
        get_cell_status_color = spec.get_cell_status_color

        # Test various warning statuses
        warning_statuses = [
            "High CV (5.2%)",
            "Low n",
            "High range (1.5)",
        ]

        for status in warning_statuses:
            color = get_cell_status_color(status)
            assert color == "#fff3cd", (
                f"Warning status '{status}' should map to yellow (#fff3cd)"
            )

    def test_get_cell_status_color_maps_errors_to_red(self, mock_streamlit):
        """
        get_cell_status_color should map error statuses to red color.

        Red (#f8d7da) indicates critical quality issues that may
        invalidate the data (e.g., Has outlier, High range > 2.0).
        """
        spec = import_module("streamlit qpcr analysis v1")
        get_cell_status_color = spec.get_cell_status_color

        # Test various error statuses
        error_statuses = [
            "Has outlier",
            "High range (2.5)",
        ]

        for status in error_statuses:
            color = get_cell_status_color(status)
            assert color == "#f8d7da", (
                f"Error status '{status}' should map to red (#f8d7da)"
            )

    def test_get_cell_display_text_formats_compact_string(self, mock_streamlit):
        """
        get_cell_display_text should format cell data into a compact
        display string like "n=3, CV=2.1%".

        This string is displayed in the grid cell and should be concise
        to fit in the cell without wrapping.
        """
        spec = import_module("streamlit qpcr analysis v1")
        get_cell_display_text = spec.get_cell_display_text

        cell_data = {
            "mean_ct": 18.5,
            "cv": 2.1,
            "status": "OK",
            "n": 3,
        }

        text = get_cell_display_text(cell_data)

        assert isinstance(text, str), "Display text must be a string"
        assert "n=" in text, "Display text should include sample count"
        assert "CV=" in text, "Display text should include CV percentage"
        assert "3" in text, "Display text should contain the n value"
        assert "2.1" in text, "Display text should contain the CV value"

    def test_get_cell_display_text_handles_edge_case_single_replicate(
        self, mock_streamlit
    ):
        """
        get_cell_display_text should handle edge case of single replicate (n=1).

        With only one replicate, CV cannot be calculated (would be 0 or undefined).
        The function should handle this gracefully.
        """
        spec = import_module("streamlit qpcr analysis v1")
        get_cell_display_text = spec.get_cell_display_text

        cell_data = {
            "mean_ct": 20.0,
            "cv": 0.0,  # CV undefined for n=1
            "status": "Low n",
            "n": 1,
        }

        text = get_cell_display_text(cell_data)

        assert isinstance(text, str), "Display text must be a string"
        assert len(text) > 0, "Display text should not be empty"
        assert "1" in text, "Display text should show n=1"

    def test_build_grid_matrix_handles_empty_dataframe(self, mock_streamlit):
        """
        build_grid_matrix should handle empty DataFrame gracefully.

        When no data is available, it should return an empty dictionary
        rather than raising an exception.
        """
        spec = import_module("streamlit qpcr analysis v1")
        build_grid_matrix = spec.build_grid_matrix
        import pandas as pd

        empty_df = pd.DataFrame()

        matrix = build_grid_matrix(empty_df)

        assert isinstance(matrix, dict), "Should return a dictionary"
        assert len(matrix) == 0, "Empty input should produce empty matrix"

    def test_build_grid_matrix_handles_single_gene_single_sample(self, mock_streamlit):
        """
        build_grid_matrix should handle minimal data: single gene, single sample.

        This edge case validates that the function works with minimal
        but valid input data.
        """
        spec = import_module("streamlit qpcr analysis v1")
        build_grid_matrix = spec.build_grid_matrix
        import pandas as pd

        # Create minimal triplicate data
        minimal_data = pd.DataFrame(
            {
                "Sample": ["Sample1", "Sample1", "Sample1"],
                "Target": ["GENE1", "GENE1", "GENE1"],
                "Mean_CT": [20.0, 20.0, 20.0],
                "SD": [0.1, 0.1, 0.1],
                "n": [3, 3, 3],
                "CV_pct": [0.5, 0.5, 0.5],
                "Range": [0.2, 0.2, 0.2],
                "Status": ["OK", "OK", "OK"],
                "Severity": ["ok", "ok", "ok"],
            }
        )

        matrix = build_grid_matrix(minimal_data)

        assert isinstance(matrix, dict), "Should return a dictionary"
        assert "GENE1" in matrix, "Should contain the gene"
        assert "Sample1" in matrix["GENE1"], "Should contain the sample"
