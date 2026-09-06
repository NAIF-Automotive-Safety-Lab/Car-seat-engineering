from stiff_medium.cone_diamond_cell_selection import (
    DIAMOND_EDGES,
    admissible_graphs,
    assess_diamond_cell_selection,
    graph_signature,
    minimal_graphs,
    powerset_edges,
)


def test_graph_scan_covers_all_four_node_edge_subsets():
    assert len(powerset_edges()) == 64


def test_diamond_is_unique_minimal_graph_under_full_constraints():
    selected = minimal_graphs(
        admissible_graphs(
            require_two_anchors=True,
            require_direct_exchange=True,
        )
    )

    assert selected == (DIAMOND_EDGES,)
    assert graph_signature(selected[0]) == ("L-A", "L-B", "L-T", "T-A", "T-B")


def test_dropping_direct_exchange_selects_square_control():
    selected = minimal_graphs(
        admissible_graphs(
            require_two_anchors=True,
            require_direct_exchange=False,
        )
    )

    assert graph_signature(selected[0]) == ("L-A", "L-B", "T-A", "T-B")
    assert len(selected[0]) == 4


def test_dropping_two_anchor_requirement_selects_one_anchor_wedge():
    selected = minimal_graphs(
        admissible_graphs(
            require_two_anchors=False,
            require_direct_exchange=True,
        )
    )
    signatures = {graph_signature(edges) for edges in selected}

    assert ("L-A", "L-T", "T-A") in signatures
    assert ("L-B", "L-T", "T-B") in signatures
    assert len(selected[0]) == 3


def test_selection_assessment_keeps_remaining_assumptions_explicit():
    result = assess_diamond_cell_selection()

    assert result.total_graphs_scanned == 64
    assert result.selected_min_edge_count == 5
    assert result.selected_min_graph_count == 1
    assert result.selected_is_diamond
    assert result.diamond_unique_under_constraints
    assert result.without_direct_min_edge_count == 4
    assert result.without_two_anchors_min_edge_count == 3
    assert not result.fully_derived
    assert "two saturated anchors plus direct branch exchange" in result.verdict
