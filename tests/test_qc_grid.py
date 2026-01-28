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
