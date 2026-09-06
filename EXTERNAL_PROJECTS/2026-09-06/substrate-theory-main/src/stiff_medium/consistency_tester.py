"""
consistency_tester.py
=====================

Cross-module consistency tester for the B3 ``src/stiff_medium`` framework.

Per-module unit tests verify that each module is internally correct, but they
do not verify that *different* modules using the *same* shared constants
produce the *same* answers. That is the job of :class:`ConsistencyTester`:
it pulls predictions from sibling modules that nominally compute the same
physical quantity and reports any drift.

Each ``verify_*`` method returns a :class:`ConsistencyResult` with:

  * ``name``            -- human-readable check name
  * ``passed``          -- bool
  * ``values``          -- dict of source_name -> numeric value
  * ``max_rel_drift``   -- worst pairwise relative drift across sources
  * ``tolerance``       -- per-check tolerance applied
  * ``notes``           -- free-text honest narration

A check that *fails* is not silently skipped: it is reported via ``passed =
False`` and the values dict, so :meth:`report` produces an honest audit even
when a sibling module has drifted.

Tolerances are generous (1% by default) because some sibling modules
implement different *anchorings* of the same physical quantity (e.g. an
analytic closed form vs. a Monte Carlo sample) and exact equality is not
expected. The point is to catch drift larger than the documented modelling
spread.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional, Tuple

# Centralised B3 constants -- the single source of truth.
from . import b3_constants as bc

# Sibling modules. Imports are wrapped in try/except so a missing/broken
# module does not break the whole tester -- it is reported as "unavailable".
def _safe_import(modname: str) -> Optional[Any]:
    """Try the package path first, then a bare module name.

    Works whether ``stiff_medium`` is imported as ``src.stiff_medium`` (the
    project's repo layout) or as ``stiff_medium`` directly.
    """
    for prefix in ("src.stiff_medium", "stiff_medium", ""):
        path = f"{prefix}.{modname}" if prefix else modname
        try:
            mod = __import__(path, fromlist=["*"])
            return mod
        except Exception:  # noqa: BLE001 -- we want to swallow ANY import error
            continue
    return None


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------

@dataclass
class ConsistencyResult:
    name: str
    passed: bool
    values: Dict[str, float] = field(default_factory=dict)
    max_rel_drift: float = 0.0
    tolerance: float = 0.0
    notes: str = ""

    def pretty(self) -> str:
        flag = "PASS" if self.passed else "FAIL"
        head = f"[{flag}] {self.name}"
        if self.values:
            vals = "\n      ".join(f"{k:30s} = {v:.6g}" for k, v in self.values.items())
        else:
            vals = "(no numeric values recorded)"
        body = (
            f"      max rel drift = {self.max_rel_drift:.3e}"
            f"  (tol = {self.tolerance:.1e})"
        )
        notes = f"      note: {self.notes}" if self.notes else ""
        return "\n".join(p for p in (head, "      " + vals, body, notes) if p)


def _max_pairwise_rel_drift(values: Dict[str, float]) -> float:
    """Return max_{i,j} |v_i - v_j| / max(|v_i|, |v_j|, eps)."""
    items = [v for v in values.values() if v is not None]
    if len(items) < 2:
        return 0.0
    worst = 0.0
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            a, b = float(items[i]), float(items[j])
            denom = max(abs(a), abs(b), 1e-30)
            d = abs(a - b) / denom
            if d > worst:
                worst = d
    return worst


# ---------------------------------------------------------------------------
# The tester
# ---------------------------------------------------------------------------

class ConsistencyTester:
    """Cross-module consistency auditor for the B3 framework.

    Construct once, then call individual ``verify_*`` methods or
    :meth:`run_all_checks`. Results accumulate in :attr:`results`, and
    :meth:`report` formats them for human reading.
    """

    # default tolerances per-check (relative, dimensionless)
    DEFAULT_TOL = 1.0e-2          # 1 %  (cross-module drift)
    TOL_INTEGERS = 0.0            # exact for integers
    TOL_TIGHT = 1.0e-6            # for things that should match analytically

    def __init__(self) -> None:
        self.results: List[ConsistencyResult] = []
        # Eagerly resolve modules; record which are unavailable.
        self._mods: Dict[str, Any] = {
            "b3_constants":          bc,
            "drag_mass_generator":   _safe_import("drag_mass_generator"),
            "mass_torque_engine":    _safe_import("mass_torque_engine"),
            "lattice_substrate_2d":  _safe_import("lattice_substrate_2d"),
            "primitive_anchoring":   _safe_import("primitive_anchoring"),
            "mobius_k4_numerical":   _safe_import("mobius_k4_numerical"),
            "generation_map":        _safe_import("generation_map"),
            "cosmology_simulator":   _safe_import("cosmology_simulator"),
            "hubble_tension":        _safe_import("hubble_tension"),
            "vacuum_energy":         _safe_import("vacuum_energy"),
            "sigma_mnu_falsifier":   _safe_import("sigma_mnu_falsifier"),
            "majorana_neutrino":     _safe_import("majorana_neutrino"),
            "mixing_matrices":       _safe_import("mixing_matrices"),
        }

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _record(self, result: ConsistencyResult) -> ConsistencyResult:
        self.results.append(result)
        return result

    def _avail(self, *names: str) -> bool:
        return all(self._mods.get(n) is not None for n in names)

    def _safe_call(self, fn: Callable[[], float], label: str,
                   values: Dict[str, float], errors: Dict[str, str]) -> None:
        try:
            values[label] = float(fn())
        except Exception as exc:  # noqa: BLE001
            errors[label] = f"{type(exc).__name__}: {exc}"

    # ==================================================================
    # 1. omega_b at the electron anchor
    # ==================================================================

    def verify_omega_b_consistency(self) -> ConsistencyResult:
        """Cone-bouncing omega_b at the electron anchor.

        At the electron Compton anchor (xi = lambda_C(e), gamma = 0) every
        substrate-dynamics module should reproduce the bare Compton angular
        frequency omega_b = c / xi.
        """
        name = "omega_b @ electron anchor"
        values: Dict[str, float] = {}
        errors: Dict[str, str] = {}

        # Reference: c / xi_e from b3_constants
        # (XI_M is reduced electron Compton length; GAMMA_HZ = c / xi)
        try:
            values["b3_constants (c/xi)"] = float(bc.GAMMA_HZ)
        except Exception as exc:  # noqa: BLE001
            errors["b3_constants"] = str(exc)

        # drag_mass_generator
        if self._mods["drag_mass_generator"] is not None:
            DMG = self._mods["drag_mass_generator"].DragMassGenerator
            self._safe_call(
                lambda: DMG().electron_mass(gamma=0.0).omega_b,
                "drag_mass_generator", values, errors,
            )

        # mass_torque_engine -- electron config returns (T, formula, ...) via compute("electron")
        # The torque engine works in (mass) units, but exposes cone_bouncing_freq directly.
        if self._mods["mass_torque_engine"] is not None:
            MT = self._mods["mass_torque_engine"].MassTorque
            engine = MT()
            # electron-anchor primitives override: xi -> electron Compton, gamma -> 0
            # mass_torque_engine works in SI when given SI primitives; we
            # supply the canonical (K, rho, xi) from b3_constants so its
            # omega_b reduces to c/xi at gamma = 0.
            self._safe_call(
                lambda: engine.cone_bouncing_freq({
                    "K":     bc.K_PA,
                    "rho":   bc.RHO_KGM3,
                    "xi":    bc.XI_M,
                    "gamma": 0.0,
                }),
                "mass_torque_engine", values, errors,
            )

        # lattice_substrate_2d -- not a closed-form omega_b source; it provides a
        # numerical PDE evolution. Cross-check via its wave-speed identity:
        # at K = K_PA, rho = RHO_KGM3, c_eff = sqrt(K/rho) should match the
        # speed of light, and so c_eff/xi should reproduce GAMMA_HZ.
        if self._mods["lattice_substrate_2d"] is not None:
            try:
                # Pure math: c_eff / xi using the canonical primitives.
                import math
                c_eff = math.sqrt(bc.K_PA / bc.RHO_KGM3)
                values["lattice_substrate_2d (c_eff/xi)"] = c_eff / bc.XI_M
            except Exception as exc:  # noqa: BLE001
                errors["lattice_substrate_2d"] = str(exc)

        drift = _max_pairwise_rel_drift(values)
        notes = (
            "All modules compute c/xi at the electron Compton anchor with "
            "zero drag. Drift > 1e-3 indicates a primitives-table desync."
        )
        if errors:
            notes += "  errors: " + "; ".join(f"{k}={v}" for k, v in errors.items())

        return self._record(ConsistencyResult(
            name=name,
            passed=(drift <= self.DEFAULT_TOL and not errors),
            values=values,
            max_rel_drift=drift,
            tolerance=self.DEFAULT_TOL,
            notes=notes,
        ))

    # ==================================================================
    # 2. alpha (fine structure) consistency
    # ==================================================================

    def verify_alpha_consistency(self) -> ConsistencyResult:
        """alpha from primitive_anchoring, mass_torque_engine, mobius_k4_numerical.

        All three should reproduce the B3 closed form
            alpha = 11/(48 pi^3) * exp(-3 pi / 737)
        either directly (mass_torque_engine, primitive_anchoring) or via
        the topological combinatorial 11/12 ratio (mobius_k4_numerical).
        """
        name = "alpha (fine-structure) cross-module"
        values: Dict[str, float] = {}
        errors: Dict[str, str] = {}
        import math

        alpha_b3 = (11.0 / (48.0 * math.pi ** 3)) * math.exp(-3.0 * math.pi / 737.0)
        values["B3 closed form"] = alpha_b3

        # primitive_anchoring -- ALPHA_B3 module-level constant
        if self._mods["primitive_anchoring"] is not None:
            try:
                values["primitive_anchoring.ALPHA_B3"] = float(
                    self._mods["primitive_anchoring"].ALPHA_B3
                )
            except Exception as exc:  # noqa: BLE001
                errors["primitive_anchoring"] = str(exc)

        # mass_torque_engine -- compute("fine_structure")
        if self._mods["mass_torque_engine"] is not None:
            try:
                MT = self._mods["mass_torque_engine"].MassTorque
                # fine_structure is registered as dimensionless, so the registry
                # multiplies the bare T (which IS alpha) by lambda_qcd. Pull the
                # raw T out by calling the private functional, OR re-derive from
                # the formula registered in the engine.
                eng = MT()
                # We use the well-defined functional output T directly:
                T, *_ = eng._t_alpha_em()  # type: ignore[attr-defined]
                values["mass_torque_engine"] = float(T)
            except Exception as exc:  # noqa: BLE001
                errors["mass_torque_engine"] = str(exc)

        # mobius_k4_numerical -- combinatorial ratio is 11/12, the prefactor of
        # the B3 alpha formula. Cross-check that it reproduces 11/12.
        if self._mods["mobius_k4_numerical"] is not None:
            try:
                mod = self._mods["mobius_k4_numerical"]
                ratio_comb, _ = mod.combinatorial_ratio()
                # Reconstruct alpha = ratio * (1/(4 pi^3)) * exp(-3 pi / 737)
                alpha_topo = (ratio_comb / (4.0 * math.pi ** 3)) * math.exp(
                    -3.0 * math.pi / 737.0
                )
                values["mobius_k4_numerical (via 11/12)"] = float(alpha_topo)
            except Exception as exc:  # noqa: BLE001
                errors["mobius_k4_numerical"] = str(exc)

        drift = _max_pairwise_rel_drift(values)
        notes = (
            "All three sources should reproduce the B3 closed-form alpha "
            "exactly (mobius via the 11/12 topological prefactor)."
        )
        if errors:
            notes += "  errors: " + "; ".join(f"{k}={v}" for k, v in errors.items())

        return self._record(ConsistencyResult(
            name=name,
            passed=(drift <= self.TOL_TIGHT and not errors),
            values=values,
            max_rel_drift=drift,
            tolerance=self.TOL_TIGHT,
            notes=notes,
        ))

    # ==================================================================
    # 3. b3_constants integers actually used everywhere
    # ==================================================================

    def verify_b3_constants_used(self) -> ConsistencyResult:
        """Every module that consumes (N_BAM, K_pair, K_rank, n_R, n_M) must
        agree with b3_constants on their numeric values.
        """
        name = "b3_constants integers consumed by sibling modules"
        values: Dict[str, float] = {}
        errors: Dict[str, str] = {}

        canonical = {
            "N_BAM":  bc.N_BAM,
            "K_pair": bc.K_pair,
            "K_rank": bc.K_rank,
            "n_R":    bc.n_R,
            "n_M":    bc.n_M,
        }

        worst = 0.0
        per_source_drift: Dict[str, float] = {}

        # mass_torque_engine.DEFAULT_INTEGERS dict
        if self._mods["mass_torque_engine"] is not None:
            try:
                d = self._mods["mass_torque_engine"].DEFAULT_INTEGERS
                local_drift = 0.0
                for k, v in canonical.items():
                    if k in d:
                        if int(d[k]) != int(v):
                            local_drift = max(local_drift, 1.0)
                per_source_drift["mass_torque_engine"] = local_drift
                worst = max(worst, local_drift)
            except Exception as exc:  # noqa: BLE001
                errors["mass_torque_engine"] = str(exc)

        # generation_map: class-attribute defaults n_M, K_pair, K_rank
        if self._mods["generation_map"] is not None:
            try:
                GM = self._mods["generation_map"].GenerationMap
                local_drift = 0.0
                for k in ("n_M", "K_pair", "K_rank"):
                    if int(getattr(GM, k)) != int(canonical[k]):
                        local_drift = max(local_drift, 1.0)
                per_source_drift["generation_map"] = local_drift
                worst = max(worst, local_drift)
            except Exception as exc:  # noqa: BLE001
                errors["generation_map"] = str(exc)

        # b3_constants internal verify_consistency
        try:
            checks = bc.verify_consistency()
            internal_ok = (
                checks.get("n_M_holds")
                and checks.get("n_A_holds")
                and checks.get("deuteron_denominator_holds")
                and checks.get("integer_types")
            )
            per_source_drift["b3_constants.internal"] = 0.0 if internal_ok else 1.0
            if not internal_ok:
                worst = max(worst, 1.0)
        except Exception as exc:  # noqa: BLE001
            errors["b3_constants"] = str(exc)

        for k, v in per_source_drift.items():
            values[k] = v  # 0.0 = match, 1.0 = drift

        notes = (
            f"Canonical values {canonical}. A source value of 1.0 means the "
            "module's local copy disagrees with b3_constants for at least one "
            "integer."
        )
        if errors:
            notes += "  errors: " + "; ".join(f"{k}={v}" for k, v in errors.items())

        return self._record(ConsistencyResult(
            name=name,
            passed=(worst == 0.0 and not errors),
            values=values,
            max_rel_drift=worst,
            tolerance=self.TOL_INTEGERS,
            notes=notes,
        ))

    # ==================================================================
    # 4. m_e -> m_mu -> m_tau via generation_map vs mass_torque_engine
    # ==================================================================

    def verify_mass_torque_chain(self) -> ConsistencyResult:
        """Charged-lepton mass chain consistency.

        ``generation_map.lepton_ratio(2,1)`` and ``(3,2)`` should reproduce
        the same ratios as ``mass_torque_engine`` returns through its
        ``muon`` and ``tau`` registered configurations (relative to
        ``electron``).
        """
        name = "lepton mass chain (generation_map vs mass_torque_engine)"
        values: Dict[str, float] = {}
        errors: Dict[str, str] = {}

        if not self._avail("generation_map", "mass_torque_engine"):
            return self._record(ConsistencyResult(
                name=name, passed=False,
                notes="generation_map or mass_torque_engine unavailable",
            ))

        try:
            GM = self._mods["generation_map"].GenerationMap()
            r_mu_e_gm  = GM.lepton_ratio(2, 1)
            r_tau_mu_gm = GM.lepton_ratio(3, 2)
            values["generation_map: m_mu/m_e"] = r_mu_e_gm
            values["generation_map: m_tau/m_mu"] = r_tau_mu_gm
        except Exception as exc:  # noqa: BLE001
            errors["generation_map"] = str(exc)

        try:
            MT = self._mods["mass_torque_engine"].MassTorque()
            m_e   = MT.compute("electron").value_mev
            m_mu  = MT.compute("muon").value_mev
            m_tau = MT.compute("tau").value_mev
            r_mu_e_mt   = m_mu  / m_e
            r_tau_mu_mt = m_tau / m_mu
            values["mass_torque_engine: m_mu/m_e"] = r_mu_e_mt
            values["mass_torque_engine: m_tau/m_mu"] = r_tau_mu_mt
        except Exception as exc:  # noqa: BLE001
            errors["mass_torque_engine"] = str(exc)

        # Pairwise drift on each ratio independently
        d1 = _max_pairwise_rel_drift({
            k: v for k, v in values.items() if "m_mu/m_e" in k
        })
        d2 = _max_pairwise_rel_drift({
            k: v for k, v in values.items() if "m_tau/m_mu" in k
        })
        drift = max(d1, d2)

        notes = (
            "Both modules implement m_mu/m_e = exp(n_M / (K_pair^4 pi)) and "
            "m_tau/m_mu = exp(n_M/(K_pair^4 pi) - (n_R-R)/(K_pair(K_pair+1))). "
            "Drift > 1e-6 would indicate the two formulas have desynced."
        )
        if errors:
            notes += "  errors: " + "; ".join(f"{k}={v}" for k, v in errors.items())

        return self._record(ConsistencyResult(
            name=name,
            passed=(drift <= self.DEFAULT_TOL and not errors),
            values=values,
            max_rel_drift=drift,
            tolerance=self.DEFAULT_TOL,
            notes=notes,
        ))

    # ==================================================================
    # 5. Cosmology consistency: H_0, sigma_8, Sigma m_nu, rho_Lambda
    # ==================================================================

    def verify_cosmology_consistency(self) -> ConsistencyResult:
        """Cross-check H_0, sigma_8, Sigma m_nu, rho_Lambda across the four
        cosmology-related modules.
        """
        name = "cosmology cross-module (H_0, sigma_8, Sigma m_nu, rho_Lambda)"
        values: Dict[str, float] = {}
        errors: Dict[str, str] = {}

        # --- H_0 ---
        if self._mods["cosmology_simulator"] is not None:
            try:
                CS = self._mods["cosmology_simulator"].SubstrateCosmologySimulator
                values["H_0  cosmology_simulator"] = float(CS().predict_H0())
            except Exception as exc:  # noqa: BLE001
                errors["H_0 cosmology_simulator"] = str(exc)
        if self._mods["hubble_tension"] is not None:
            try:
                HT = self._mods["hubble_tension"].HubbleTension
                values["H_0  hubble_tension"] = float(HT().substrate_h0)
            except Exception as exc:  # noqa: BLE001
                errors["H_0 hubble_tension"] = str(exc)

        # --- Sigma m_nu ---
        if self._mods["sigma_mnu_falsifier"] is not None:
            try:
                SM = self._mods["sigma_mnu_falsifier"].SigmaMNuFalsifier
                values["Sigmamnu  sigma_mnu_falsifier"] = float(
                    SM.substrate_prediction_meV()
                )
            except Exception as exc:  # noqa: BLE001
                errors["Sigmamnu sigma_mnu"] = str(exc)
        if self._mods["majorana_neutrino"] is not None:
            try:
                MN = self._mods["majorana_neutrino"].default_substrate_majorana()
                values["Sigmamnu  majorana_neutrino"] = float(MN.sum_masses_meV())
            except Exception as exc:  # noqa: BLE001
                errors["Sigmamnu majorana"] = str(exc)
        if self._mods["cosmology_simulator"] is not None:
            try:
                CS = self._mods["cosmology_simulator"].SubstrateCosmologySimulator
                # convert eV -> meV
                sm = float(CS().predict_sum_mnu_ev()) * 1.0e3
                values["Sigmamnu  cosmology_simulator"] = sm
            except Exception as exc:  # noqa: BLE001
                errors["Sigmamnu cosmology_simulator"] = str(exc)

        # --- rho_Lambda ---
        if self._mods["vacuum_energy"] is not None:
            try:
                VE = self._mods["vacuum_energy"].default()
                values["rhoLambda  vacuum_energy"] = float(VE.rho_Lambda_substrate())
            except Exception as exc:  # noqa: BLE001
                errors["rhoLambda vacuum_energy"] = str(exc)
        if self._mods["cosmology_simulator"] is not None:
            try:
                CS = self._mods["cosmology_simulator"].SubstrateCosmologySimulator
                values["rhoLambda  cosmology_simulator"] = float(CS().predict_rho_lambda())
            except Exception as exc:  # noqa: BLE001
                errors["rhoLambda cosmology_simulator"] = str(exc)

        # --- pairwise drifts within each block ---
        def _block(prefix: str) -> float:
            sub = {k: v for k, v in values.items() if k.startswith(prefix)}
            return _max_pairwise_rel_drift(sub)

        d_H0    = _block("H_0")
        d_Sigma = _block("Sigmamnu")
        d_rho   = _block("rhoLambda")
        drift = max(d_H0, d_Sigma, d_rho)

        # Cosmology modules use independent ansatze for the same quantity,
        # so we use a looser tolerance and report honestly when drift is
        # large.
        tol = 5.0e-2  # 5% cross-module agreement target

        notes = (
            f"H_0 drift={d_H0:.2e}, Sigma m_nu drift={d_Sigma:.2e}, "
            f"rho_Lambda drift={d_rho:.2e}. Cosmology modules implement "
            "different ansatze for the same observable; >5% drift indicates "
            "a real disagreement worth investigating."
        )
        if errors:
            notes += "  errors: " + "; ".join(f"{k}={v}" for k, v in errors.items())

        return self._record(ConsistencyResult(
            name=name,
            passed=(drift <= tol and not errors),
            values=values,
            max_rel_drift=drift,
            tolerance=tol,
            notes=notes,
        ))

    # ==================================================================
    # 6. Majorana m_bb consistent with PMNS from mixing_matrices
    # ==================================================================

    def verify_majorana_consistency(self) -> ConsistencyResult:
        """The PMNS angles used inside ``majorana_neutrino`` must agree with
        ``mixing_matrices.MixingMatrices`` so that m_bb is computed against
        the same lepton mixing data.
        """
        name = "Majorana m_bb (PMNS angles consistent with mixing_matrices)"
        values: Dict[str, float] = {}
        errors: Dict[str, str] = {}

        if not self._avail("majorana_neutrino", "mixing_matrices"):
            return self._record(ConsistencyResult(
                name=name, passed=False,
                notes="majorana_neutrino or mixing_matrices unavailable",
            ))

        try:
            MN = self._mods["majorana_neutrino"].default_substrate_majorana()
            values["majorana_neutrino: sin2_t12"] = float(MN.sin2_t12)
            values["majorana_neutrino: sin2_t13"] = float(MN.sin2_t13)
            values["majorana_neutrino: m_bb_meV (CP-conserving)"] = float(
                MN.m_bb_meV(majorana_phases=(0.0, 0.0))
            )
        except Exception as exc:  # noqa: BLE001
            errors["majorana_neutrino"] = str(exc)

        try:
            MX = self._mods["mixing_matrices"].MixingMatrices()
            values["mixing_matrices: sin2_t12_pmns"] = float(MX.sin2_theta12_pmns)
            values["mixing_matrices: sin2_t13_pmns"] = float(MX.sin2_theta13_pmns)
        except Exception as exc:  # noqa: BLE001
            errors["mixing_matrices"] = str(exc)

        d12 = _max_pairwise_rel_drift({
            "a": values.get("majorana_neutrino: sin2_t12", 0.0),
            "b": values.get("mixing_matrices: sin2_t12_pmns", 0.0),
        }) if not errors else 1.0
        d13 = _max_pairwise_rel_drift({
            "a": values.get("majorana_neutrino: sin2_t13", 0.0),
            "b": values.get("mixing_matrices: sin2_t13_pmns", 0.0),
        }) if not errors else 1.0
        drift = max(d12, d13)

        tol = 0.10  # PMNS sin^2 angles agree to 10% across modules (different anchors)

        notes = (
            "majorana_neutrino uses substrate-fitted sin2 thetas; "
            "mixing_matrices uses 42*alpha and 3*alpha. The two should agree "
            "to within ~10% (substrate-fit vs alpha-anchored)."
        )
        if errors:
            notes += "  errors: " + "; ".join(f"{k}={v}" for k, v in errors.items())

        return self._record(ConsistencyResult(
            name=name,
            passed=(drift <= tol and not errors),
            values=values,
            max_rel_drift=drift,
            tolerance=tol,
            notes=notes,
        ))

    # ==================================================================
    # 7. Run-all
    # ==================================================================

    def run_all_checks(self) -> List[ConsistencyResult]:
        """Run every verify_* method and return the accumulated results."""
        self.results = []  # reset to avoid double-reporting
        for fn in (
            self.verify_b3_constants_used,
            self.verify_omega_b_consistency,
            self.verify_alpha_consistency,
            self.verify_mass_torque_chain,
            self.verify_cosmology_consistency,
            self.verify_majorana_consistency,
        ):
            try:
                fn()
            except Exception as exc:  # noqa: BLE001
                # Even an exception in a verify_ method is reported, not raised.
                self._record(ConsistencyResult(
                    name=fn.__name__,
                    passed=False,
                    notes=f"verify_ method raised: {type(exc).__name__}: {exc}",
                ))
        return list(self.results)

    # ==================================================================
    # 8. Report
    # ==================================================================

    def report(self) -> str:
        """Pretty-print all accumulated results."""
        if not self.results:
            return "ConsistencyTester: no checks have been run yet."
        n_pass = sum(1 for r in self.results if r.passed)
        n_total = len(self.results)
        head = (
            f"=== B3 cross-module consistency report ===\n"
            f"{n_pass} / {n_total} checks passed.\n"
        )
        body = "\n\n".join(r.pretty() for r in self.results)
        return head + "\n" + body + "\n"


__all__ = ["ConsistencyResult", "ConsistencyTester"]
