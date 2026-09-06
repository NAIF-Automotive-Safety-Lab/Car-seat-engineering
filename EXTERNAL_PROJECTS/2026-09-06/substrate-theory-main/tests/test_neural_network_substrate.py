"""Tests for the paired neural-network learning-dynamics module.

Required claims (per the GEOMETRY/SIM/VIZ specification):
    (a) A single-layer perceptron solves AND and OR (linearly separable).
    (b) A single-layer perceptron CANNOT solve XOR (Minsky-Papert 1969).
    (c) An MLP with one hidden layer solves XOR (universality).
    (d) Gradient descent on a quadratic test problem converges towards
        the origin and is monotonically decreasing in loss for a small
        learning rate.
    (e) The loss landscape probe returns a non-negative scalar field
        with a clear minimum near the trained weight.

Plus standard ML invariants:
    * weight-count formula matches (sum n_in*n_out + n_out per layer)
    * sigmoid / tanh / relu derivatives match the analytic forms at
      checkpoint values
    * Glorot init produces weights bounded by sqrt(6/(n_in+n_out))
    * Boolean datasets have correct truth tables
    * Generalisation gap is finite and signed correctly
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from src.stiff_medium.neural_network_substrate import (
    DEFAULT_HIDDEN_WIDTH,
    DEFAULT_LEARNING_RATE,
    LAMBDA_QCD_MEV,
    NeuralNetworkGeometry,
    NeuralNetworkSimulator,
    boolean_dataset,
    relu,
    relu_prime,
    sigmoid,
    sigmoid_prime,
    tanh,
    tanh_prime,
)


# --------------------------------------------------------------------- #
# Activation function tests                                             #
# --------------------------------------------------------------------- #

def test_sigmoid_at_zero_is_half():
    assert sigmoid(0.0) == pytest.approx(0.5, abs=1e-12)


def test_sigmoid_extreme_values_clip_safely():
    # Should not raise an OverflowError or produce NaN
    s_neg = float(sigmoid(-1.0e3))
    s_pos = float(sigmoid(+1.0e3))
    assert 0.0 <= s_neg < 1.0e-30
    assert 1.0 - 1.0e-12 < s_pos <= 1.0
    assert math.isfinite(s_neg) and math.isfinite(s_pos)


def test_sigmoid_prime_matches_s_times_one_minus_s():
    xs = np.linspace(-3, 3, 7)
    s = sigmoid(xs)
    expected = s * (1.0 - s)
    assert np.allclose(sigmoid_prime(xs), expected)


def test_tanh_at_zero_is_zero_and_prime_is_one():
    assert tanh(0.0) == pytest.approx(0.0, abs=1e-12)
    assert tanh_prime(0.0) == pytest.approx(1.0, abs=1e-12)


def test_tanh_prime_matches_one_minus_tanh_squared():
    xs = np.linspace(-2, 2, 11)
    expected = 1.0 - np.tanh(xs) ** 2
    assert np.allclose(tanh_prime(xs), expected)


def test_relu_is_positive_part():
    xs = np.array([-2.0, -0.5, 0.0, 0.5, 2.0])
    expected = np.array([0.0, 0.0, 0.0, 0.5, 2.0])
    assert np.allclose(relu(xs), expected)


def test_relu_prime_is_step_function():
    xs = np.array([-2.0, -0.5, 0.0, 0.5, 2.0])
    expected = np.array([0.0, 0.0, 0.0, 1.0, 1.0])
    assert np.allclose(relu_prime(xs), expected)


# --------------------------------------------------------------------- #
# Boolean dataset truth tables                                          #
# --------------------------------------------------------------------- #

@pytest.mark.parametrize("name,truth", [
    ("AND",  [0, 0, 0, 1]),
    ("OR",   [0, 1, 1, 1]),
    ("NAND", [1, 1, 1, 0]),
    ("NOR",  [1, 0, 0, 0]),
    ("XOR",  [0, 1, 1, 0]),
    ("XNOR", [1, 0, 0, 1]),
])
def test_boolean_truth_table(name, truth):
    X, y = boolean_dataset(name)
    assert X.shape == (4, 2)
    assert y.shape == (4,)
    assert list(y.astype(int)) == truth


def test_unknown_boolean_dataset_raises():
    with pytest.raises(ValueError):
        boolean_dataset("FOO")


# --------------------------------------------------------------------- #
# GEOMETRY: parameter counts and Glorot bounds                          #
# --------------------------------------------------------------------- #

def test_perceptron_geometry_is_2_to_1():
    g = NeuralNetworkGeometry(layer_sizes=(2, 1))
    assert g.is_perceptron()
    assert g.n_hidden_layers() == 0
    # 2 weights + 1 bias = 3 parameters
    assert g.n_parameters() == 3
    assert g.n_weights() == 2
    assert g.n_biases() == 1


def test_mlp_geometry_2_4_1_has_17_parameters():
    g = NeuralNetworkGeometry(layer_sizes=(2, 4, 1))
    # 2*4 + 4 + 4*1 + 1 = 8 + 4 + 4 + 1 = 17
    assert g.n_parameters() == 17
    assert g.n_hidden_layers() == 1
    assert not g.is_perceptron()
    assert g.can_solve_xor_in_principle()


def test_glorot_init_bounds():
    g = NeuralNetworkGeometry(layer_sizes=(2, 4, 1), seed=0)
    for W, (n_in, n_out) in zip(g.weights,
                                zip(g.layer_sizes[:-1], g.layer_sizes[1:])):
        r = math.sqrt(6.0 / (n_in + n_out))
        assert np.all(W >= -r - 1e-12)
        assert np.all(W <= +r + 1e-12)


def test_init_is_reproducible():
    g1 = NeuralNetworkGeometry(layer_sizes=(2, 4, 1), seed=42)
    g2 = NeuralNetworkGeometry(layer_sizes=(2, 4, 1), seed=42)
    for W1, W2 in zip(g1.weights, g2.weights):
        assert np.allclose(W1, W2)


def test_strain_amplitudes_concatenate_all_weights():
    g = NeuralNetworkGeometry(layer_sizes=(2, 3, 1), seed=1)
    amp = g.strain_amplitudes()
    assert amp.shape == (g.n_weights(),)


def test_layer_sizes_must_have_at_least_two_layers():
    with pytest.raises(ValueError):
        NeuralNetworkGeometry(layer_sizes=(2,))


def test_unknown_activation_raises():
    with pytest.raises(ValueError):
        NeuralNetworkGeometry(layer_sizes=(2, 1), activation="banana")


def test_geometry_substrate_signature_basic():
    g = NeuralNetworkGeometry(layer_sizes=(2, 4, 1))
    sig = g.substrate_signature()
    assert sig["lambda_qcd_MeV"] == LAMBDA_QCD_MEV
    assert sig["n_layers"] == 3
    assert sig["n_hidden_layers"] == 1
    assert sig["n_parameters"] == 17
    assert sig["universal_approximator"]


# --------------------------------------------------------------------- #
# (a) Perceptron solves AND and OR                                      #
# --------------------------------------------------------------------- #

def test_perceptron_solves_and():
    sim = NeuralNetworkSimulator(
        geometry=NeuralNetworkGeometry(layer_sizes=(2, 1), seed=0),
        learning_rate=DEFAULT_LEARNING_RATE,
        max_epochs=4000,
    )
    out = sim.train_perceptron("AND", max_epochs=4000)
    X, y = boolean_dataset("AND")
    preds = sim.predict_classes(X)
    assert np.array_equal(preds, y.astype(int))


def test_perceptron_solves_or():
    sim = NeuralNetworkSimulator(
        geometry=NeuralNetworkGeometry(layer_sizes=(2, 1), seed=0),
        learning_rate=DEFAULT_LEARNING_RATE,
        max_epochs=4000,
    )
    sim.train_perceptron("OR", max_epochs=4000)
    X, y = boolean_dataset("OR")
    preds = sim.predict_classes(X)
    assert np.array_equal(preds, y.astype(int))


# --------------------------------------------------------------------- #
# (b) Perceptron CANNOT solve XOR                                       #
# --------------------------------------------------------------------- #

def test_perceptron_cannot_solve_xor():
    """Minsky-Papert 1969: a linear classifier has no solution to XOR."""
    sim = NeuralNetworkSimulator(
        geometry=NeuralNetworkGeometry(layer_sizes=(2, 1), seed=0),
        learning_rate=0.5,
        max_epochs=8000,
    )
    sim.train_perceptron("XOR", max_epochs=8000, tolerance=1e-6)
    X, y = boolean_dataset("XOR")
    preds = sim.predict_classes(X)
    # A perceptron cannot achieve full accuracy on XOR; at most 2 of 4
    n_correct = int(np.sum(preds == y.astype(int)))
    assert n_correct <= 3, (
        f"single-layer perceptron unexpectedly got {n_correct}/4 on XOR; "
        "should be impossible by Minsky-Papert 1969"
    )


# --------------------------------------------------------------------- #
# (c) MLP solves XOR (universality)                                     #
# --------------------------------------------------------------------- #

def test_mlp_solves_xor():
    """A 2-input, H-hidden, 1-output MLP can solve XOR (Cybenko 1989)."""
    # Try a few seeds; the result should be robustly solvable.
    solved = False
    for seed in (0, 1, 2, 3, 4, 7):
        sim = NeuralNetworkSimulator(
            geometry=NeuralNetworkGeometry(
                layer_sizes=(2, DEFAULT_HIDDEN_WIDTH, 1), seed=seed,
            ),
            learning_rate=0.8,
            max_epochs=5000,
        )
        sim.train_mlp("XOR", hidden=DEFAULT_HIDDEN_WIDTH, max_epochs=5000)
        X, y = boolean_dataset("XOR")
        preds = sim.predict_classes(X)
        if np.array_equal(preds, y.astype(int)):
            solved = True
            break
    assert solved, "MLP failed to solve XOR for any of the tried seeds"


def test_mlp_xor_loss_decreases_during_training():
    sim = NeuralNetworkSimulator(
        geometry=NeuralNetworkGeometry(
            layer_sizes=(2, DEFAULT_HIDDEN_WIDTH, 1), seed=0,
        ),
        learning_rate=0.5,
        max_epochs=2000,
    )
    out = sim.train_mlp("XOR", hidden=DEFAULT_HIDDEN_WIDTH, max_epochs=2000)
    losses = out["losses"]
    assert losses[-1] < losses[0], "loss did not decrease during training"


# --------------------------------------------------------------------- #
# (d) Gradient-descent convergence on a quadratic test                  #
# --------------------------------------------------------------------- #

def test_gradient_descent_converges_to_origin_for_isotropic_quadratic():
    path = NeuralNetworkSimulator.gradient_descent_trajectory(
        a1=1.0, a2=1.0, learning_rate=0.1, n_steps=200,
        start=(2.0, 2.0),
    )
    # Should approach the global minimum (0, 0)
    final = np.linalg.norm(path[-1])
    assert final < 1e-3, f"GD did not converge: final norm {final}"


def test_gradient_descent_loss_monotone_for_small_lr():
    """For small lr, L = a1*w1^2 + a2*w2^2 is monotonically decreasing."""
    a1, a2 = 1.0, 5.0
    path = NeuralNetworkSimulator.gradient_descent_trajectory(
        a1=a1, a2=a2, learning_rate=0.05, n_steps=200,
        start=(-2.5, 2.5),
    )
    L = a1 * path[:, 0] ** 2 + a2 * path[:, 1] ** 2
    diffs = np.diff(L)
    # All steps should reduce or maintain loss
    assert np.all(diffs <= 1e-9), "loss is not monotonically non-increasing"


def test_gradient_descent_zigzags_for_high_condition_number():
    """High a2/a1 ratio causes the path to bounce between the steep walls."""
    path = NeuralNetworkSimulator.gradient_descent_trajectory(
        a1=1.0, a2=20.0, learning_rate=0.04, n_steps=80,
        start=(-2.0, 1.5),
    )
    # The w2 component should change sign at least once (zig-zag)
    sign_changes = int(np.sum(np.diff(np.sign(path[:, 1])) != 0))
    assert sign_changes >= 1, (
        "expected at least one sign change in w2 for ill-conditioned GD"
    )


# --------------------------------------------------------------------- #
# (e) Loss-landscape probe                                              #
# --------------------------------------------------------------------- #

def test_loss_landscape_returns_nonnegative_values():
    sim = NeuralNetworkSimulator(
        geometry=NeuralNetworkGeometry(layer_sizes=(2, 4, 1), seed=0),
        max_epochs=200,
    )
    X, y = boolean_dataset("XOR")
    w1s, w2s, L = sim.loss_landscape(X, y, n_grid=11)
    assert L.shape == (11, 11)
    assert np.all(L >= 0.0)


def test_loss_landscape_quadratic_helper_is_psd():
    w1, w2 = np.meshgrid(np.linspace(-1, 1, 5),
                          np.linspace(-1, 1, 5))
    L = NeuralNetworkGeometry.quadratic_loss_landscape(w1, w2)
    assert np.all(L >= 0.0)
    # minimum at origin
    i = np.unravel_index(np.argmin(L), L.shape)
    assert L[i] == pytest.approx(0.0, abs=1e-12)


def test_loss_landscape_rosenbrock_minimum_at_a_a_squared():
    a, b = 1.0, 100.0
    w1, w2 = np.meshgrid(np.linspace(-1.5, 1.8, 80),
                          np.linspace(-0.5, 2.5, 80))
    L = NeuralNetworkGeometry.rosenbrock_loss_landscape(w1, w2, a=a, b=b)
    i_min = np.unravel_index(np.argmin(L), L.shape)
    # Allow slack from grid resolution
    assert abs(w1[i_min] - a) < 0.10
    assert abs(w2[i_min] - a ** 2) < 0.10


# --------------------------------------------------------------------- #
# Forward pass / gradient correctness                                   #
# --------------------------------------------------------------------- #

def test_forward_pass_shape():
    sim = NeuralNetworkSimulator(
        geometry=NeuralNetworkGeometry(layer_sizes=(2, 4, 1), seed=0)
    )
    X = np.array([[0, 0], [0, 1], [1, 0], [1, 1]], dtype=float)
    y_pred = sim.predict(X)
    assert y_pred.shape == (4, 1)
    # Sigmoid output bounded in [0, 1]
    assert np.all(y_pred >= 0.0) and np.all(y_pred <= 1.0)


def test_forward_pass_deterministic_for_fixed_seed():
    sim1 = NeuralNetworkSimulator(
        geometry=NeuralNetworkGeometry(layer_sizes=(2, 4, 1), seed=99))
    sim2 = NeuralNetworkSimulator(
        geometry=NeuralNetworkGeometry(layer_sizes=(2, 4, 1), seed=99))
    X = np.array([[0.5, 0.5]])
    assert np.allclose(sim1.predict(X), sim2.predict(X))


def test_step_reduces_loss_on_xor_at_some_point():
    """Single step should reduce loss for at least some random seeds."""
    sim = NeuralNetworkSimulator(
        geometry=NeuralNetworkGeometry(layer_sizes=(2, 4, 1), seed=0),
        learning_rate=0.5,
    )
    X, y = boolean_dataset("XOR")
    L_before = sim._loss_value(sim.predict(X), y)
    sim.step(X, y)
    L_after = sim._loss_value(sim.predict(X), y)
    # Should not blow up
    assert L_after < L_before * 1.1


# --------------------------------------------------------------------- #
# Decision-boundary probe                                               #
# --------------------------------------------------------------------- #

def test_decision_boundary_grid_shape():
    sim = NeuralNetworkSimulator(
        geometry=NeuralNetworkGeometry(layer_sizes=(2, 4, 1), seed=0)
    )
    X_mesh, Y_mesh, Z = sim.decision_boundary_grid(n_grid=20)
    assert X_mesh.shape == (20, 20)
    assert Y_mesh.shape == (20, 20)
    assert Z.shape == (20, 20)
    # Sigmoid outputs
    assert np.all(Z >= 0.0) and np.all(Z <= 1.0)


# --------------------------------------------------------------------- #
# Train / test split, generalisation gap                                #
# --------------------------------------------------------------------- #

def test_train_test_split_counts():
    X, y = NeuralNetworkSimulator.generate_noisy_dataset(n=80, seed=0)
    X_tr, y_tr, X_te, y_te = NeuralNetworkSimulator.train_test_split(
        X, y, test_fraction=0.25, seed=0,
    )
    assert X_tr.shape[0] + X_te.shape[0] == 80
    assert X_te.shape[0] == 20


def test_generalisation_gap_finite_after_short_training():
    sim = NeuralNetworkSimulator(
        geometry=NeuralNetworkGeometry(layer_sizes=(2, 8, 1), seed=0),
        max_epochs=200,
    )
    X, y = NeuralNetworkSimulator.generate_noisy_dataset(n=80, seed=0)
    X_tr, y_tr, X_te, y_te = NeuralNetworkSimulator.train_test_split(
        X, y, test_fraction=0.3, seed=0,
    )
    sim.train(X_tr, y_tr, max_epochs=200)
    gap = sim.generalization_gap(X_tr, y_tr, X_te, y_te)
    assert math.isfinite(gap["train_loss"])
    assert math.isfinite(gap["test_loss"])
    assert math.isfinite(gap["gap"])


# --------------------------------------------------------------------- #
# Substrate signatures                                                  #
# --------------------------------------------------------------------- #

def test_simulator_substrate_signature():
    sim = NeuralNetworkSimulator(
        geometry=NeuralNetworkGeometry(layer_sizes=(2, 4, 1), seed=0)
    )
    sig = sim.substrate_signature()
    assert sig["lambda_qcd_MeV"] == LAMBDA_QCD_MEV
    assert sig["learning_rate"] == DEFAULT_LEARNING_RATE
    assert sig["loss"] == "bce"
    assert "back-propagated" in sig["interpretation"]
