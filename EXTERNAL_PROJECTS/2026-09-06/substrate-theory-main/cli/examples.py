"""Ten worked examples of the substrate_physics library.

Run directly:

    python examples.py

Each example prints the value, the measured value (when known), the
precision estimate, and the substrate module where the formula lives.
"""

from __future__ import annotations

from substrate_physics import SubstratePhysics, format_prediction


def banner(idx: int, title: str) -> None:
    print()
    print(f"=== Example {idx}: {title} ===")


def main() -> None:
    sp = SubstratePhysics()

    banner(1, "Fine-structure constant α from substrate (closed form)")
    p = sp.predict_alpha_em()
    print(format_prediction(p))
    print(f"  inv_alpha_substrate = {1.0 / p.value:.6f}")

    banner(2, "Muon/electron mass ratio from inventory integers")
    p = sp.predict_lepton_mass_ratio()
    print(format_prediction(p))

    banner(3, "Tau/electron mass ratio from lepton tower")
    p = sp.predict_tau_over_electron_ratio()
    print(format_prediction(p))

    banner(4, "Atomic ionisation energies (H..Ar, light elements)")
    for sym in ("H", "C", "N", "O", "Si"):
        print(format_prediction(sp.predict_atomic_ie(sym)))

    banner(5, "Semiconductor bandgaps")
    for mat in ("silicon", "germanium", "diamond", "gaas"):
        print(format_prediction(sp.predict_bandgap(mat)))

    banner(6, "Hadron masses from face-spin v4")
    for h in ("p", "n", "Lambda", "Delta", "Omega", "JPsi"):
        print(format_prediction(sp.predict_hadron_mass(h)))

    banner(7, "Lifetimes (V-A weak-coupling inheritance)")
    for part in ("muon", "tau", "neutron", "pi+"):
        print(format_prediction(sp.predict_lifetime(part)))

    banner(8, "Cosmology / BSM predictions")
    print(format_prediction(sp.predict_dark_matter_mass()))
    print(format_prediction(sp.predict_neutrino_sum()))
    print(format_prediction(sp.predict_hierarchy()))

    banner(9, "Materials science: T_c ceiling, Madelung, Debye, BCS")
    print(format_prediction(sp.predict_tc_max()))
    print(format_prediction(sp.predict_madelung("nacl")))
    print(format_prediction(sp.predict_debye_temperature("diamond")))
    print(format_prediction(sp.predict_bcs_gap_ratio("weak_coupling")))
    print(format_prediction(sp.predict_grueneisen("copper")))

    banner(10, "Engineering: fracture plastic-zone radius")
    # K_I = 50 MPa·sqrt(m), sigma_y = 600 MPa  -> typical mild-steel sample
    print(format_prediction(sp.predict_fracture_zone(K_I=50.0, sigma_y=600.0)))


if __name__ == "__main__":
    main()
