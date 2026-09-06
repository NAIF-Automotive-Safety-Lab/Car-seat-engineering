"""API documentation generator for the stiff_medium framework.

Introspection-driven: scans the package directory, imports modules
(best-effort), and extracts class / function signatures, docstrings,
and inter-module imports. Produces markdown reference docs and a
graphviz dependency graph.

Usage:
    from stiff_medium.api_documentation import APIDocumentation
    api = APIDocumentation()
    api.report()
    api.markdown_export("API.md")
    api.dependency_graph_dot("deps.dot")
"""

from __future__ import annotations

import ast
import importlib
import inspect
import os
import pkgutil
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

PACKAGE_NAME = "stiff_medium"

# Modules considered "entry points" — high-level orchestrators / main APIs.
ENTRY_POINT_HINTS = (
    "MassTorque",
    "GenerationMap",
    "BaryonFaceSpin",
    "CalibrationEngine",
    "ConsistencyTester",
    "CosmologySimulator",
    "DependencyLedger",
    "B3Constants",
    "AlphaDerivation",
    "BlackHole",
)

# Classes that have shown up frequently in usage / integration tests
# — used by usage_examples() to seed template snippets.
TOP_USED_CLASSES = (
    "MassTorque",
    "GenerationMap",
    "BaryonFaceSpin",
    "CalibrationEngine",
    "ConsistencyTester",
    "CosmologySimulator",
    "AlphaDerivation",
    "BlackHole",
    "B3Constants",
    "DependencyLedger",
)


@dataclass
class SymbolInfo:
    name: str
    kind: str  # "class" | "function"
    signature: str
    docstring: str
    module: str


@dataclass
class ModuleInfo:
    name: str  # short name (no package prefix)
    qualname: str  # full dotted path
    path: str
    docstring: str = ""
    symbols: List[SymbolInfo] = field(default_factory=list)
    imports: List[str] = field(default_factory=list)
    importable: bool = True
    import_error: str = ""


class APIDocumentation:
    """Introspection-based API documentation generator."""

    def __init__(self, package: str = PACKAGE_NAME, package_dir: Optional[str] = None):
        self.package = package
        if package_dir is None:
            try:
                pkg = importlib.import_module(package)
                package_dir = os.path.dirname(pkg.__file__)
            except Exception:
                # Fall back to assuming src/ layout next to cwd.
                here = Path(__file__).resolve().parent
                package_dir = str(here)
        self.package_dir = package_dir
        self._modules: Optional[List[ModuleInfo]] = None

    # ------------------------------------------------------------------
    # 1. list_all_modules
    # ------------------------------------------------------------------
    def list_all_modules(self) -> List[str]:
        """Return sorted list of public module short names in the package."""
        names: List[str] = []
        for entry in sorted(os.listdir(self.package_dir)):
            if not entry.endswith(".py"):
                continue
            if entry.startswith("_"):
                continue
            names.append(entry[:-3])
        return names

    # ------------------------------------------------------------------
    # internal: collect ModuleInfo for everything (cached)
    # ------------------------------------------------------------------
    def _collect(self) -> List[ModuleInfo]:
        if self._modules is not None:
            return self._modules
        infos: List[ModuleInfo] = []
        for short in self.list_all_modules():
            qual = f"{self.package}.{short}"
            path = os.path.join(self.package_dir, short + ".py")
            info = ModuleInfo(name=short, qualname=qual, path=path)
            # AST-based docstring + imports (always works, no execution).
            try:
                with open(path, "r", encoding="utf-8") as fh:
                    source = fh.read()
                tree = ast.parse(source)
                info.docstring = ast.get_docstring(tree) or ""
                info.imports = self._extract_imports(tree)
            except Exception as exc:  # pragma: no cover - defensive
                info.docstring = ""
                info.imports = []
                info.importable = False
                info.import_error = f"parse: {exc!r}"
            infos.append(info)
        self._modules = infos
        return infos

    @staticmethod
    def _extract_imports(tree: ast.AST) -> List[str]:
        out: List[str] = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    out.append(alias.name)
            elif isinstance(node, ast.ImportFrom):
                if node.module:
                    out.append(node.module)
        return out

    # ------------------------------------------------------------------
    # 2. module_signatures
    # ------------------------------------------------------------------
    def module_signatures(self, module_name: str) -> List[Dict[str, str]]:
        """Extract class/function signatures + docstrings from a module.

        Tries live import first (gives true signatures). Falls back to AST
        if the module fails to import (common for modules with optional
        heavy deps).
        """
        qual = (
            module_name
            if module_name.startswith(self.package + ".")
            else f"{self.package}.{module_name}"
        )
        results: List[Dict[str, str]] = []
        try:
            mod = importlib.import_module(qual)
        except Exception:
            mod = None

        if mod is not None:
            for name, obj in inspect.getmembers(mod):
                if name.startswith("_"):
                    continue
                if not (inspect.isclass(obj) or inspect.isfunction(obj)):
                    continue
                # Restrict to objects defined in this module (skip re-exports).
                obj_mod = getattr(obj, "__module__", "")
                if obj_mod != qual:
                    continue
                try:
                    sig = str(inspect.signature(obj))
                except (TypeError, ValueError):
                    sig = "(...)"
                doc = inspect.getdoc(obj) or ""
                kind = "class" if inspect.isclass(obj) else "function"
                results.append(
                    {
                        "name": name,
                        "kind": kind,
                        "signature": f"{name}{sig}",
                        "docstring": doc,
                        "module": qual,
                    }
                )
            return results

        # AST fallback.
        path = os.path.join(self.package_dir, module_name.split(".")[-1] + ".py")
        try:
            with open(path, "r", encoding="utf-8") as fh:
                tree = ast.parse(fh.read())
        except Exception:
            return results
        for node in tree.body:
            if isinstance(node, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
                if node.name.startswith("_"):
                    continue
                kind = "class" if isinstance(node, ast.ClassDef) else "function"
                args = self._format_args(node) if not isinstance(node, ast.ClassDef) else "(...)"
                results.append(
                    {
                        "name": node.name,
                        "kind": kind,
                        "signature": f"{node.name}{args}",
                        "docstring": ast.get_docstring(node) or "",
                        "module": qual,
                    }
                )
        return results

    @staticmethod
    def _format_args(node: ast.AST) -> str:
        try:
            return "(" + ", ".join(a.arg for a in node.args.args) + ")"
        except Exception:
            return "(...)"

    # ------------------------------------------------------------------
    # 3. cross_module_imports
    # ------------------------------------------------------------------
    def cross_module_imports(self) -> Dict[str, List[str]]:
        """Return {module_name: [internal modules it imports]} graph."""
        infos = self._collect()
        known = {info.name for info in infos}
        graph: Dict[str, List[str]] = {}
        for info in infos:
            internal: List[str] = []
            for imp in info.imports:
                # Strip package prefix if present.
                tail = imp.split(".")
                # Look for a token matching a known internal module.
                for tok in tail:
                    if tok in known and tok != info.name:
                        internal.append(tok)
                        break
                else:
                    # Also handle "stiff_medium.foo" → foo
                    if tail and tail[0] == self.package and len(tail) > 1:
                        cand = tail[1]
                        if cand in known and cand != info.name:
                            internal.append(cand)
            graph[info.name] = sorted(set(internal))
        return graph

    # ------------------------------------------------------------------
    # 4. entry_points
    # ------------------------------------------------------------------
    def entry_points(self) -> List[Dict[str, str]]:
        """Identify modules whose top-level classes look like main entry points."""
        out: List[Dict[str, str]] = []
        for info in self._collect():
            try:
                with open(info.path, "r", encoding="utf-8") as fh:
                    tree = ast.parse(fh.read())
            except Exception:
                continue
            for node in tree.body:
                if isinstance(node, ast.ClassDef) and node.name in ENTRY_POINT_HINTS:
                    out.append(
                        {
                            "module": info.name,
                            "class": node.name,
                            "docstring": (ast.get_docstring(node) or "").splitlines()[:1][0]
                            if ast.get_docstring(node)
                            else "",
                        }
                    )
        return out

    # ------------------------------------------------------------------
    # 5. api_completeness
    # ------------------------------------------------------------------
    def api_completeness(self) -> Dict[str, Any]:
        """Classify modules by whether they expose a public API."""
        infos = self._collect()
        public, internal, undocumented = [], [], []
        for info in infos:
            try:
                with open(info.path, "r", encoding="utf-8") as fh:
                    tree = ast.parse(fh.read())
            except Exception:
                internal.append(info.name)
                continue
            public_syms = [
                n
                for n in tree.body
                if isinstance(n, (ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef))
                and not n.name.startswith("_")
            ]
            if public_syms:
                public.append(info.name)
            else:
                internal.append(info.name)
            if not ast.get_docstring(tree):
                undocumented.append(info.name)
        return {
            "total": len(infos),
            "public_count": len(public),
            "internal_count": len(internal),
            "undocumented_count": len(undocumented),
            "public": sorted(public),
            "internal": sorted(internal),
            "undocumented": sorted(undocumented),
            "coverage_pct": (
                100.0 * len(public) / len(infos) if infos else 0.0
            ),
        }

    # ------------------------------------------------------------------
    # 6. usage_examples
    # ------------------------------------------------------------------
    def usage_examples(self) -> Dict[str, str]:
        """Return small code snippets for the top-used classes."""
        examples: Dict[str, str] = {}
        # Build {class_name: module_name} index from AST.
        index: Dict[str, str] = {}
        for info in self._collect():
            try:
                with open(info.path, "r", encoding="utf-8") as fh:
                    tree = ast.parse(fh.read())
            except Exception:
                continue
            for node in tree.body:
                if isinstance(node, ast.ClassDef) and not node.name.startswith("_"):
                    index.setdefault(node.name, info.name)

        for cls in TOP_USED_CLASSES[:10]:
            if cls not in index:
                continue
            mod = index[cls]
            examples[cls] = (
                f"from {self.package}.{mod} import {cls}\n"
                f"obj = {cls}()\n"
                f"# Inspect available methods:\n"
                f"# print([m for m in dir(obj) if not m.startswith('_')])\n"
            )
        return examples

    # ------------------------------------------------------------------
    # 7. markdown_export
    # ------------------------------------------------------------------
    def markdown_export(self, filename: str) -> str:
        """Write API docs to `filename` as markdown. Returns the path."""
        infos = self._collect()
        completeness = self.api_completeness()
        entries = self.entry_points()
        graph = self.cross_module_imports()
        examples = self.usage_examples()

        lines: List[str] = []
        lines.append(f"# {self.package} — API Reference\n")
        lines.append(
            f"_Auto-generated by `APIDocumentation`. "
            f"{len(infos)} modules, {completeness['public_count']} with public API "
            f"({completeness['coverage_pct']:.1f}% coverage)._\n"
        )

        lines.append("\n## Table of Contents\n")
        lines.append("- [Overview](#overview)")
        lines.append("- [Entry Points](#entry-points)")
        lines.append("- [Usage Examples](#usage-examples)")
        lines.append("- [Modules](#modules)")
        lines.append("- [Dependency Graph](#dependency-graph)\n")

        lines.append("\n## Overview\n")
        lines.append(f"- Total modules: **{completeness['total']}**")
        lines.append(f"- Public-API modules: **{completeness['public_count']}**")
        lines.append(f"- Internal-only modules: **{completeness['internal_count']}**")
        lines.append(
            f"- Undocumented modules: **{completeness['undocumented_count']}**\n"
        )

        lines.append("\n## Entry Points\n")
        if entries:
            for e in entries:
                lines.append(
                    f"- `{e['class']}` (in `{self.package}.{e['module']}`) — {e['docstring']}"
                )
        else:
            lines.append("_No entry-point classes detected._")

        lines.append("\n\n## Usage Examples\n")
        for cls, snippet in examples.items():
            lines.append(f"### {cls}\n")
            lines.append("```python")
            lines.append(snippet.rstrip())
            lines.append("```\n")

        lines.append("\n## Modules\n")
        for info in infos:
            lines.append(f"\n### `{info.qualname}`\n")
            if info.docstring:
                lines.append(info.docstring.strip().splitlines()[0])
                lines.append("")
            sigs = self.module_signatures(info.name)
            if not sigs:
                lines.append("_(no public symbols)_\n")
                continue
            for s in sigs:
                lines.append(f"- **{s['kind']}** `{s['signature']}`")
                if s["docstring"]:
                    first = s["docstring"].strip().splitlines()[0]
                    lines.append(f"  - {first}")
            deps = graph.get(info.name, [])
            if deps:
                lines.append(f"\n_Internal deps: {', '.join('`' + d + '`' for d in deps)}_")

        lines.append("\n\n## Dependency Graph\n")
        lines.append(
            "See accompanying `.dot` file (generated by `dependency_graph_dot`)."
        )

        out = "\n".join(lines) + "\n"
        with open(filename, "w", encoding="utf-8") as fh:
            fh.write(out)
        return filename

    # ------------------------------------------------------------------
    # 8. dependency_graph_dot
    # ------------------------------------------------------------------
    def dependency_graph_dot(self, filename: str) -> str:
        """Write a graphviz .dot file of internal module dependencies."""
        graph = self.cross_module_imports()
        entries = {e["module"] for e in self.entry_points()}
        # Also include nodes for the major class-bearing modules.
        major_class_nodes = set()
        for info in self._collect():
            try:
                with open(info.path, "r", encoding="utf-8") as fh:
                    tree = ast.parse(fh.read())
            except Exception:
                continue
            for node in tree.body:
                if isinstance(node, ast.ClassDef) and node.name in ENTRY_POINT_HINTS:
                    major_class_nodes.add(info.name)

        lines: List[str] = []
        lines.append("digraph stiff_medium {")
        lines.append('  rankdir="LR";')
        lines.append('  node [shape=box, fontname="Helvetica", fontsize=10];')
        for mod in sorted(graph):
            attrs = []
            if mod in entries or mod in major_class_nodes:
                attrs.append('style=filled')
                attrs.append('fillcolor="#ffd966"')
            attrstr = (" [" + ", ".join(attrs) + "]") if attrs else ""
            lines.append(f'  "{mod}"{attrstr};')
        for mod, deps in sorted(graph.items()):
            for d in deps:
                lines.append(f'  "{mod}" -> "{d}";')
        lines.append("}")
        out = "\n".join(lines) + "\n"
        with open(filename, "w", encoding="utf-8") as fh:
            fh.write(out)
        return filename

    # ------------------------------------------------------------------
    # 9. report
    # ------------------------------------------------------------------
    def report(self) -> str:
        """Human-readable summary string. Also prints to stdout."""
        infos = self._collect()
        comp = self.api_completeness()
        entries = self.entry_points()
        graph = self.cross_module_imports()
        # Most-imported (in-degree) modules.
        in_deg: Counter = Counter()
        for deps in graph.values():
            for d in deps:
                in_deg[d] += 1
        most_used = in_deg.most_common(10)

        lines = []
        lines.append("=" * 64)
        lines.append(f" {self.package}: API documentation report")
        lines.append("=" * 64)
        lines.append(f"Modules scanned        : {comp['total']}")
        lines.append(f"With public API        : {comp['public_count']}  ({comp['coverage_pct']:.1f}%)")
        lines.append(f"Internal only          : {comp['internal_count']}")
        lines.append(f"Undocumented (no top-  : {comp['undocumented_count']}")
        lines.append(f"  level docstring)")
        lines.append("")
        lines.append(f"Entry-point classes detected: {len(entries)}")
        for e in entries:
            lines.append(f"  - {e['class']:<24s} in {self.package}.{e['module']}")
        lines.append("")
        lines.append("Top in-degree (most-imported) internal modules:")
        for name, n in most_used:
            lines.append(f"  - {name:<32s} {n} importers")
        lines.append("=" * 64)
        text = "\n".join(lines)
        print(text)
        return text


if __name__ == "__main__":  # pragma: no cover
    APIDocumentation().report()
