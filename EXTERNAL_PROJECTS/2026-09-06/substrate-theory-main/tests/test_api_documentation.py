"""Tests for stiff_medium.api_documentation.APIDocumentation."""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

# Ensure src/ is importable when tests are run from repo root.
ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from stiff_medium.api_documentation import APIDocumentation  # noqa: E402


@pytest.fixture(scope="module")
def api() -> APIDocumentation:
    return APIDocumentation()


def test_list_all_modules_has_at_least_40(api: APIDocumentation) -> None:
    mods = api.list_all_modules()
    assert isinstance(mods, list)
    assert all(isinstance(m, str) for m in mods)
    assert len(mods) >= 40, f"expected >=40 modules, got {len(mods)}"
    # No private modules.
    assert all(not m.startswith("_") for m in mods)
    # Sorted.
    assert mods == sorted(mods)


def test_list_modules_contains_known_examples(api: APIDocumentation) -> None:
    mods = set(api.list_all_modules())
    # A handful of modules we know exist in this codebase.
    expected_any = {
        "b3_constants",
        "calibration_engine",
        "consistency_tester",
        "cosmology_simulator",
    }
    assert expected_any & mods, f"expected at least one of {expected_any} in modules"


def test_module_signatures_returns_list(api: APIDocumentation) -> None:
    mods = api.list_all_modules()
    # Find any module that yields signatures (some may be empty/internal).
    saw_signatures = False
    for m in mods[:30]:
        sigs = api.module_signatures(m)
        assert isinstance(sigs, list)
        if sigs:
            saw_signatures = True
            s = sigs[0]
            assert {"name", "kind", "signature", "docstring", "module"} <= set(s.keys())
            assert s["kind"] in {"class", "function"}
    assert saw_signatures, "no module yielded any public signatures (unexpected)"


def test_cross_module_imports_is_graph(api: APIDocumentation) -> None:
    g = api.cross_module_imports()
    assert isinstance(g, dict)
    assert len(g) >= 40
    for k, v in g.items():
        assert isinstance(v, list)
        assert k not in v  # no self-loops


def test_entry_points_detected(api: APIDocumentation) -> None:
    eps = api.entry_points()
    assert isinstance(eps, list)
    # We expect at least *some* entry-point class to exist in this codebase.
    # Don't hard-require a specific one — just structural sanity.
    for e in eps:
        assert {"module", "class", "docstring"} <= set(e.keys())


def test_api_completeness_shape(api: APIDocumentation) -> None:
    c = api.api_completeness()
    assert c["total"] >= 40
    assert c["public_count"] + c["internal_count"] == c["total"]
    assert 0.0 <= c["coverage_pct"] <= 100.0
    assert isinstance(c["public"], list)
    assert isinstance(c["internal"], list)


def test_usage_examples_returns_dict(api: APIDocumentation) -> None:
    ex = api.usage_examples()
    assert isinstance(ex, dict)
    assert len(ex) <= 10
    for cls, snippet in ex.items():
        assert isinstance(cls, str)
        assert "import" in snippet
        assert cls in snippet


def test_markdown_export_creates_file_with_sections(
    api: APIDocumentation, tmp_path: Path
) -> None:
    out = tmp_path / "API.md"
    path = api.markdown_export(str(out))
    assert path == str(out)
    assert out.exists()
    text = out.read_text()
    # Required sections.
    for section in (
        "# stiff_medium",
        "## Overview",
        "## Entry Points",
        "## Usage Examples",
        "## Modules",
    ):
        assert section in text, f"missing section: {section!r}"
    # Has substantial content (40+ modules → file should be sizeable).
    assert len(text) > 5000


def test_dependency_graph_dot_has_major_class_nodes(
    api: APIDocumentation, tmp_path: Path
) -> None:
    out = tmp_path / "deps.dot"
    api.dependency_graph_dot(str(out))
    assert out.exists()
    text = out.read_text()
    assert text.startswith("digraph")
    assert text.rstrip().endswith("}")
    # Should declare nodes.
    assert '"' in text
    # Should reference at least one major class-bearing module.
    mods = set(api.list_all_modules())
    major_candidates = {
        "b3_constants",
        "calibration_engine",
        "consistency_tester",
        "cosmology_simulator",
        "dependency_ledger",
    } & mods
    assert any(f'"{m}"' in text for m in major_candidates), (
        "dependency graph missing nodes for major modules"
    )
    # Should contain at least one edge (-> arrow), given 40+ interconnected modules.
    assert "->" in text


def test_report_returns_string(api: APIDocumentation, capsys) -> None:
    text = api.report()
    assert isinstance(text, str)
    assert "API documentation report" in text
    captured = capsys.readouterr()
    assert "API documentation report" in captured.out
