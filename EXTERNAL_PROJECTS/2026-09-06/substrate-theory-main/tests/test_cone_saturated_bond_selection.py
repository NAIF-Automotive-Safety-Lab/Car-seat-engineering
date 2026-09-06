from stiff_medium.cone_saturated_bond_selection import (
    assess_saturated_bond_selection,
    critical_core_cost_for_single_bond,
    distributed_slip_energy,
    saturation_barrier_energy,
    selected_bond_count,
    slip_energy_table,
)


def test_saturation_barrier_energy_is_convex_near_cap():
    assert saturation_barrier_energy(0.0) == 0.0
    assert saturation_barrier_energy(0.9) > saturation_barrier_energy(0.5)
    assert saturation_barrier_energy(-0.9) == saturation_barrier_energy(0.9)


def test_pure_barrier_delocalizes_fixed_phase_slip():
    table = slip_energy_table(0.9, (1, 2, 4, 8, 16, 32, 64))
    energies = [energy for _, energy in table]

    assert selected_bond_count(table) == 64
    assert energies[0] > energies[-1]


def test_positive_core_cost_can_select_one_bond_but_must_be_derived():
    counts = (1, 2, 4, 8, 16, 32, 64)
    critical = critical_core_cost_for_single_bond(0.9, counts)
    table = slip_energy_table(0.9, counts, core_cost_per_active_bond=1.01 * critical)

    assert critical > 1.0
    assert selected_bond_count(table) == 1


def test_saturated_bond_selection_reports_honest_failure():
    result = assess_saturated_bond_selection()

    assert result.total_strain_fraction == 0.9
    assert result.pure_selected_bond_count == 64
    assert result.pure_one_bond_energy > result.pure_widest_energy
    assert not result.pure_barrier_selects_single_bond
    assert result.critical_core_cost > 1.0
    assert result.core_selected_bond_count == 1
    assert not result.core_cost_fixed_by_substrate
    assert not result.fully_derived
    assert "barrier alone delocalizes" in result.verdict


def test_invalid_saturated_bond_selection_inputs_are_rejected():
    for f in (1.0, -1.0):
        try:
            saturation_barrier_energy(f)
        except ValueError as exc:
            assert "strain_fraction" in str(exc)
        else:
            raise AssertionError("invalid strain fraction should fail")

    try:
        distributed_slip_energy(0.5, 0)
    except ValueError as exc:
        assert "bond_count" in str(exc)
    else:
        raise AssertionError("invalid bond count should fail")

    try:
        distributed_slip_energy(0.5, 1, core_cost_per_active_bond=-1.0)
    except ValueError as exc:
        assert "core_cost" in str(exc)
    else:
        raise AssertionError("invalid core cost should fail")

    try:
        selected_bond_count(())
    except ValueError as exc:
        assert "non-empty" in str(exc)
    else:
        raise AssertionError("empty energy table should fail")

    try:
        critical_core_cost_for_single_bond(0.5, (2, 4))
    except ValueError as exc:
        assert "include 1" in str(exc)
    else:
        raise AssertionError("missing one-bond control should fail")
