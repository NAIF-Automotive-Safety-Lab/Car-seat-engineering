"""Paired neural-network learning-dynamics GEOMETRY + SIMULATOR.

Connects the B3 / stiff-medium ontology to artificial neural networks
(perceptron, MLP, gradient descent) by interpreting the trainable weight
matrix as a discrete substrate strain pattern that is iteratively reshaped
by error back-propagation.

GEOMETRY
--------
``NeuralNetworkGeometry`` represents a feed-forward network as a layered
graph of substrate nodes connected by oriented "strain bonds" (weights).
Each weight is a local strain-gradient amplitude; the bias is a baseline
saturation offset.  The full weight tensor is the strain pattern of the
network.  The loss landscape L(theta) is then the substrate's curvature
energy as a function of weight configuration.

Standard textbook reference:
    * Rosenblatt 1958     -- the perceptron (linear separator)
    * Minsky-Papert 1969  -- single-layer perceptrons cannot solve XOR
    * Rumelhart-Hinton-
      Williams 1986       -- back-propagation in multilayer perceptrons
    * universality (Cybenko 1989; Hornik 1991): a single hidden layer of
      sufficient width with a sigmoidal activation can approximate any
      continuous function on a compact domain.

SIMULATOR
---------
``NeuralNetworkSimulator`` exposes:

    * ``train_perceptron`` -- batch gradient descent on a single linear
      perceptron with sigmoid output and binary cross-entropy loss.
    * ``train_mlp`` -- back-propagation training of a 2-input, H-hidden,
      1-output MLP with tanh hidden activation and sigmoid output.  This
      is the smallest network capable of solving XOR.
    * ``loss_landscape`` -- evaluate L(w_1, w_2) over a 2D slice of the
      MLP's weight space, used to visualise the loss surface.
    * ``gradient_descent_trajectory`` -- record the (w_1, w_2) path taken
      by SGD on a 2-parameter quadratic test problem
      L(w) = 0.5 * w^T A w (Hessian eigenvalue ratio = condition number).
    * ``train_test_split`` and ``generalization_gap`` -- toy demonstration
      of overfitting (train loss << test loss when capacity > samples).

Substrate interpretation:
    * weight w_ij = oriented strain bond between substrate nodes i, j
    * activation = local saturation response to incoming strain
    * loss L(theta) = total substrate distortion / output-mismatch energy
    * gradient descent = strain-pattern relaxation under back-propagated
      error currents (analogous to a viscous gradient-flow on the loss
      landscape, with the learning rate playing the role of inverse
      viscosity)
    * universality (single hidden layer is enough) reflects the
      completeness of substrate strain modes: any smooth target field can
      be written as a finite superposition of localised strain bumps, the
      activation functions provide the basis.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, List, Optional, Tuple

import math
import numpy as np


# ---------------------------------------------------------------------------
# Substrate / framework constants
# ---------------------------------------------------------------------------

LAMBDA_QCD_MEV: float = 200.0

# Default training hyperparameters
DEFAULT_LEARNING_RATE: float = 0.5
DEFAULT_MAX_EPOCHS: int = 2000
DEFAULT_TOLERANCE: float = 1.0e-3

# Default MLP hidden width (XOR needs >= 2; we use 4 for headroom)
DEFAULT_HIDDEN_WIDTH: int = 4


# ---------------------------------------------------------------------------
# Activation functions and their derivatives
# ---------------------------------------------------------------------------

def sigmoid(x: np.ndarray | float) -> np.ndarray:
    """Stable logistic sigmoid."""
    x = np.asarray(x, dtype=float)
    out = np.where(x >= 0,
                   1.0 / (1.0 + np.exp(-np.clip(x, -500.0, 500.0))),
                   np.exp(np.clip(x, -500.0, 500.0))
                   / (1.0 + np.exp(np.clip(x, -500.0, 500.0))))
    return out


def sigmoid_prime(x: np.ndarray | float) -> np.ndarray:
    s = sigmoid(x)
    return s * (1.0 - s)


def tanh(x: np.ndarray | float) -> np.ndarray:
    return np.tanh(x)


def tanh_prime(x: np.ndarray | float) -> np.ndarray:
    return 1.0 - np.tanh(x) ** 2


def relu(x: np.ndarray | float) -> np.ndarray:
    return np.maximum(np.asarray(x, dtype=float), 0.0)


def relu_prime(x: np.ndarray | float) -> np.ndarray:
    return (np.asarray(x, dtype=float) > 0.0).astype(float)


# ---------------------------------------------------------------------------
# Canonical Boolean datasets
# ---------------------------------------------------------------------------

def boolean_dataset(name: str) -> Tuple[np.ndarray, np.ndarray]:
    """Return (X, y) for AND, OR, NAND, NOR, XOR, XNOR.

    X has shape (4, 2) covering all 2-bit inputs; y has shape (4,).
    """
    X = np.array([[0, 0], [0, 1], [1, 0], [1, 1]], dtype=float)
    table = {
        "AND":  np.array([0, 0, 0, 1], dtype=float),
        "OR":   np.array([0, 1, 1, 1], dtype=float),
        "NAND": np.array([1, 1, 1, 0], dtype=float),
        "NOR":  np.array([1, 0, 0, 0], dtype=float),
        "XOR":  np.array([0, 1, 1, 0], dtype=float),
        "XNOR": np.array([1, 0, 0, 1], dtype=float),
    }
    if name not in table:
        raise ValueError(f"unknown dataset {name!r}; "
                         f"choose from {sorted(table)}")
    return X, table[name]


# ---------------------------------------------------------------------------
# GEOMETRY
# ---------------------------------------------------------------------------

@dataclass
class NeuralNetworkGeometry:
    """Layered substrate strain pattern: weights as oriented strain bonds.

    Attributes
    ----------
    layer_sizes : tuple[int, ...]
        Number of nodes in each layer, including input and output layers.
        For a perceptron with 2 inputs and 1 output: (2, 1).
        For an XOR-solving MLP with 4 hidden units: (2, 4, 1).
    activation : str
        Hidden-layer activation: "tanh" (default) or "sigmoid" or "relu".
    output_activation : str
        Output-layer activation: "sigmoid" (default, for binary outputs).
    seed : int
        RNG seed for reproducible weight initialisation.
    """

    layer_sizes: Tuple[int, ...] = (2, 1)
    activation: str = "tanh"
    output_activation: str = "sigmoid"
    seed: int = 0

    weights: List[np.ndarray] = field(default_factory=list)
    biases:  List[np.ndarray] = field(default_factory=list)

    def __post_init__(self) -> None:
        if len(self.layer_sizes) < 2:
            raise ValueError("layer_sizes must have at least input + output")
        if self.activation not in ("tanh", "sigmoid", "relu"):
            raise ValueError(f"activation {self.activation!r} not supported")
        if self.output_activation not in ("sigmoid", "tanh", "linear"):
            raise ValueError(
                f"output_activation {self.output_activation!r} not supported"
            )
        if not self.weights:
            self.initialize_weights()

    # ------------------------------------------------------------------
    # Initialisation
    # ------------------------------------------------------------------

    def initialize_weights(self, seed: Optional[int] = None) -> None:
        """Glorot-style initialisation: w ~ U(-r, +r), r = sqrt(6/(n_in+n_out))."""
        rng = np.random.default_rng(self.seed if seed is None else seed)
        self.weights = []
        self.biases  = []
        for n_in, n_out in zip(self.layer_sizes[:-1], self.layer_sizes[1:]):
            r = math.sqrt(6.0 / (n_in + n_out))
            W = rng.uniform(-r, +r, size=(n_in, n_out))
            b = np.zeros(n_out)
            self.weights.append(W)
            self.biases.append(b)

    # ------------------------------------------------------------------
    # Counts
    # ------------------------------------------------------------------

    def n_layers(self) -> int:
        """Total number of layers (input + hidden + output)."""
        return len(self.layer_sizes)

    def n_hidden_layers(self) -> int:
        return max(0, len(self.layer_sizes) - 2)

    def n_parameters(self) -> int:
        """Total trainable weights + biases."""
        n = 0
        for n_in, n_out in zip(self.layer_sizes[:-1], self.layer_sizes[1:]):
            n += n_in * n_out + n_out
        return n

    def n_weights(self) -> int:
        return sum(W.size for W in self.weights)

    def n_biases(self) -> int:
        return sum(b.size for b in self.biases)

    # ------------------------------------------------------------------
    # Strain-pattern view
    # ------------------------------------------------------------------

    def strain_amplitudes(self) -> np.ndarray:
        """Return the flat vector of all weights (substrate strain pattern)."""
        return np.concatenate([W.flatten() for W in self.weights])

    def strain_pattern_norm(self) -> float:
        """Frobenius norm of the flat strain pattern (weight magnitude scale)."""
        amp = self.strain_amplitudes()
        return float(np.sqrt(np.sum(amp ** 2)))

    def is_perceptron(self) -> bool:
        """True if no hidden layers (single-layer perceptron)."""
        return self.n_hidden_layers() == 0

    def can_solve_xor_in_principle(self) -> bool:
        """Universality: an MLP needs at least one hidden unit (>=2 helps)."""
        return self.n_hidden_layers() >= 1 and self.layer_sizes[1] >= 2

    # ------------------------------------------------------------------
    # Loss landscape: scalar L(theta) projected onto a 2D slice
    # ------------------------------------------------------------------

    @staticmethod
    def quadratic_loss_landscape(
        w1: np.ndarray, w2: np.ndarray,
        a1: float = 1.0, a2: float = 5.0,
        center: Tuple[float, float] = (0.0, 0.0),
    ) -> np.ndarray:
        """A 2D quadratic loss bowl L(w1, w2) = a1*(w1-c1)^2 + a2*(w2-c2)^2.

        The eigenvalue ratio a2/a1 sets the condition number; a high ratio
        produces an elongated valley that gradient descent zigzags through.
        """
        c1, c2 = center
        return a1 * (w1 - c1) ** 2 + a2 * (w2 - c2) ** 2

    @staticmethod
    def rosenbrock_loss_landscape(
        w1: np.ndarray, w2: np.ndarray, a: float = 1.0, b: float = 100.0,
    ) -> np.ndarray:
        """Rosenbrock 'banana' valley: L = (a - w1)^2 + b*(w2 - w1^2)^2.

        Classic non-quadratic test surface with a curved valley and a
        unique minimum at (w1, w2) = (a, a^2).
        """
        return (a - w1) ** 2 + b * (w2 - w1 ** 2) ** 2

    # ------------------------------------------------------------------
    # Substrate metadata
    # ------------------------------------------------------------------

    def substrate_signature(self) -> dict:
        return {
            "interpretation": (
                "neural network = layered substrate; weight = oriented "
                "strain bond; activation = local saturation response; "
                "loss surface = substrate distortion energy; gradient "
                "descent = strain-pattern relaxation flow"
            ),
            "lambda_qcd_MeV": LAMBDA_QCD_MEV,
            "layer_sizes": self.layer_sizes,
            "n_layers": self.n_layers(),
            "n_hidden_layers": self.n_hidden_layers(),
            "n_parameters": self.n_parameters(),
            "is_perceptron": self.is_perceptron(),
            "universal_approximator": self.can_solve_xor_in_principle(),
        }


# ---------------------------------------------------------------------------
# SIMULATOR
# ---------------------------------------------------------------------------

@dataclass
class NeuralNetworkSimulator:
    """Forward/backward propagation, gradient descent, loss-landscape probes.

    Trains the network in ``geometry`` in-place by mutating
    ``geometry.weights`` and ``geometry.biases``.  All training is plain
    full-batch (or stochastic) gradient descent on binary cross-entropy
    or mean-squared-error loss.
    """

    geometry: NeuralNetworkGeometry = field(default_factory=NeuralNetworkGeometry)
    learning_rate: float = DEFAULT_LEARNING_RATE
    max_epochs: int = DEFAULT_MAX_EPOCHS
    tolerance: float = DEFAULT_TOLERANCE
    loss: str = "bce"   # binary cross-entropy ("bce") or mean-squared error ("mse")

    # ------------------------------------------------------------------
    # Activation lookup
    # ------------------------------------------------------------------

    def _act(self, kind: str) -> Tuple[Callable, Callable]:
        if kind == "tanh":
            return tanh, tanh_prime
        if kind == "sigmoid":
            return sigmoid, sigmoid_prime
        if kind == "relu":
            return relu, relu_prime
        if kind == "linear":
            return (lambda x: x, lambda x: np.ones_like(x))
        raise ValueError(f"unknown activation {kind!r}")

    # ------------------------------------------------------------------
    # Forward pass
    # ------------------------------------------------------------------

    def forward(self, X: np.ndarray) -> Tuple[np.ndarray, dict]:
        """Forward propagate inputs X (shape: n_samples, n_features).

        Returns (output, cache) where cache stores pre-activations and
        post-activations for use in the backward pass.
        """
        g = self.geometry
        a = X
        zs: List[np.ndarray] = []
        as_: List[np.ndarray] = [X]
        n_layers_w = len(g.weights)
        for i, (W, b) in enumerate(zip(g.weights, g.biases)):
            z = a @ W + b
            zs.append(z)
            kind = g.output_activation if i == n_layers_w - 1 else g.activation
            act, _ = self._act(kind)
            a = act(z)
            as_.append(a)
        return a, {"zs": zs, "as": as_}

    def predict(self, X: np.ndarray) -> np.ndarray:
        """Return network output for X (shape: n_samples, n_outputs)."""
        out, _ = self.forward(X)
        return out

    def predict_classes(self, X: np.ndarray, threshold: float = 0.5) -> np.ndarray:
        """Hard 0/1 classification at the given decision threshold."""
        out = self.predict(X)
        if out.ndim > 1 and out.shape[1] == 1:
            out = out[:, 0]
        return (out >= threshold).astype(int)

    # ------------------------------------------------------------------
    # Loss functions
    # ------------------------------------------------------------------

    def _loss_value(self, y_pred: np.ndarray, y_true: np.ndarray) -> float:
        y_true = np.asarray(y_true, dtype=float).reshape(y_pred.shape)
        if self.loss == "mse":
            return float(np.mean((y_pred - y_true) ** 2))
        # binary cross-entropy with epsilon clamping
        eps = 1.0e-12
        yp = np.clip(y_pred, eps, 1.0 - eps)
        return float(-np.mean(y_true * np.log(yp) + (1.0 - y_true) * np.log(1.0 - yp)))

    def _output_grad(self, y_pred: np.ndarray, y_true: np.ndarray) -> np.ndarray:
        """Gradient of loss w.r.t. y_pred."""
        y_true = np.asarray(y_true, dtype=float).reshape(y_pred.shape)
        n = y_pred.shape[0]
        if self.loss == "mse":
            return 2.0 * (y_pred - y_true) / n
        # bce derivative
        eps = 1.0e-12
        yp = np.clip(y_pred, eps, 1.0 - eps)
        return (yp - y_true) / (yp * (1.0 - yp) * n)

    # ------------------------------------------------------------------
    # Backward pass
    # ------------------------------------------------------------------

    def backward(self, cache: dict, y_true: np.ndarray
                 ) -> Tuple[List[np.ndarray], List[np.ndarray]]:
        """Standard back-propagation; returns per-layer (dW, db)."""
        g = self.geometry
        zs   = cache["zs"]
        acts = cache["as"]
        n_layers_w = len(g.weights)

        grads_W: List[np.ndarray] = [np.zeros_like(W) for W in g.weights]
        grads_b: List[np.ndarray] = [np.zeros_like(b) for b in g.biases]

        # Output layer
        y_pred = acts[-1]
        out_kind = g.output_activation
        _, dact_out = self._act(out_kind)
        dL_dy = self._output_grad(y_pred, y_true)
        # If sigmoid + bce, the output-layer "delta" simplifies; we keep
        # the explicit chain rule for clarity and generality.
        delta = dL_dy * dact_out(zs[-1])

        for i in range(n_layers_w - 1, -1, -1):
            a_in = acts[i]                    # shape: (n, n_in)
            grads_W[i] = a_in.T @ delta       # (n_in, n_out)
            grads_b[i] = np.sum(delta, axis=0)
            if i > 0:
                kind = g.activation if i - 1 < n_layers_w - 1 else g.output_activation
                _, dact = self._act(kind)
                delta = (delta @ g.weights[i].T) * dact(zs[i - 1])
        return grads_W, grads_b

    # ------------------------------------------------------------------
    # Single GD step
    # ------------------------------------------------------------------

    def step(self, X: np.ndarray, y: np.ndarray) -> float:
        """One full-batch gradient-descent step; returns the post-step loss."""
        y_pred, cache = self.forward(X)
        loss_before = self._loss_value(y_pred, y)
        grads_W, grads_b = self.backward(cache, y)
        for i in range(len(self.geometry.weights)):
            self.geometry.weights[i] -= self.learning_rate * grads_W[i]
            self.geometry.biases[i]  -= self.learning_rate * grads_b[i]
        return loss_before

    # ------------------------------------------------------------------
    # Full training loop
    # ------------------------------------------------------------------

    def train(self, X: np.ndarray, y: np.ndarray,
              max_epochs: Optional[int] = None,
              tolerance: Optional[float] = None,
              record_every: int = 1) -> dict:
        """Full-batch gradient descent until loss < tolerance or max_epochs.

        Returns
        -------
        dict with keys:
            losses : list[float]   (one per recorded epoch)
            epochs : list[int]
            final_loss : float
            converged : bool
            n_epochs : int
        """
        n_eps = max_epochs if max_epochs is not None else self.max_epochs
        tol   = tolerance  if tolerance  is not None else self.tolerance

        losses: List[float] = []
        epochs: List[int]   = []
        L = float("inf")
        epoch = 0
        for epoch in range(n_eps):
            L = self.step(X, y)
            if epoch % record_every == 0:
                losses.append(L)
                epochs.append(epoch)
            if L < tol:
                # record final loss
                if losses[-1] != L:
                    losses.append(L)
                    epochs.append(epoch)
                break
        return {
            "losses": losses,
            "epochs": epochs,
            "final_loss": float(L),
            "converged": bool(L < tol),
            "n_epochs": epoch + 1,
        }

    # ------------------------------------------------------------------
    # Convenience: train_perceptron / train_mlp factory entry points
    # ------------------------------------------------------------------

    def train_perceptron(self, dataset: str = "AND",
                         max_epochs: Optional[int] = None,
                         tolerance: Optional[float] = None) -> dict:
        """Train a single-layer perceptron on the given Boolean dataset.

        Re-initialises the network as (2, 1) so the simulator becomes a
        true single-layer perceptron, then runs gradient descent.
        """
        self.geometry = NeuralNetworkGeometry(
            layer_sizes=(2, 1),
            activation="tanh",
            output_activation="sigmoid",
            seed=self.geometry.seed,
        )
        X, y = boolean_dataset(dataset)
        return self.train(X, y, max_epochs=max_epochs, tolerance=tolerance)

    def train_mlp(self, dataset: str = "XOR",
                  hidden: int = DEFAULT_HIDDEN_WIDTH,
                  max_epochs: Optional[int] = None,
                  tolerance: Optional[float] = None) -> dict:
        """Train a 2-input H-hidden 1-output MLP on the dataset (default XOR)."""
        self.geometry = NeuralNetworkGeometry(
            layer_sizes=(2, hidden, 1),
            activation="tanh",
            output_activation="sigmoid",
            seed=self.geometry.seed,
        )
        X, y = boolean_dataset(dataset)
        return self.train(X, y, max_epochs=max_epochs, tolerance=tolerance)

    # ------------------------------------------------------------------
    # Loss landscape probing
    # ------------------------------------------------------------------

    def loss_landscape(self, X: np.ndarray, y: np.ndarray,
                       w1_range: Tuple[float, float] = (-3.0, 3.0),
                       w2_range: Tuple[float, float] = (-3.0, 3.0),
                       n_grid: int = 41,
                       i_layer: int = 0,
                       i_unit:  int = 0,
                       j_unit:  int = 0,
                       k_unit:  int = 1) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Evaluate L(w_a, w_b) on a 2D slice of the network's weight space.

        Sweeps two scalar entries of W[i_layer], holding all other
        parameters fixed at their current values.  Useful for visualising
        the loss surface around the trained / random initial point.
        """
        g = self.geometry
        if i_layer < 0 or i_layer >= len(g.weights):
            raise ValueError("i_layer out of range")
        W_orig = g.weights[i_layer].copy()
        n_in, n_out = W_orig.shape
        if not (0 <= i_unit < n_in and 0 <= j_unit < n_out):
            raise ValueError("i_unit/j_unit out of range")
        if not (0 <= i_unit < n_in and 0 <= k_unit < n_out):
            raise ValueError("k_unit out of range")

        w1s = np.linspace(*w1_range, n_grid)
        w2s = np.linspace(*w2_range, n_grid)
        L_grid = np.zeros((n_grid, n_grid))
        try:
            for ii, w1 in enumerate(w1s):
                for jj, w2 in enumerate(w2s):
                    g.weights[i_layer][i_unit, j_unit] = w1
                    g.weights[i_layer][i_unit, k_unit] = w2
                    y_pred, _ = self.forward(X)
                    L_grid[ii, jj] = self._loss_value(y_pred, y)
        finally:
            g.weights[i_layer][:] = W_orig
        return w1s, w2s, L_grid

    # ------------------------------------------------------------------
    # Gradient-descent trajectory on a 2-parameter quadratic test
    # ------------------------------------------------------------------

    @staticmethod
    def gradient_descent_trajectory(
        a1: float = 1.0, a2: float = 5.0,
        learning_rate: float = 0.1,
        n_steps: int = 60,
        start: Tuple[float, float] = (-2.5, 2.5),
    ) -> np.ndarray:
        """Trajectory of vanilla GD on L(w) = a1*w1^2 + a2*w2^2.

        Returns array of shape (n_steps + 1, 2) with the (w1, w2) path.
        Eigenvalue ratio a2/a1 controls zig-zag intensity.
        """
        path = np.zeros((n_steps + 1, 2))
        path[0] = start
        for i in range(n_steps):
            w1, w2 = path[i]
            g1 = 2.0 * a1 * w1
            g2 = 2.0 * a2 * w2
            path[i + 1, 0] = w1 - learning_rate * g1
            path[i + 1, 1] = w2 - learning_rate * g2
        return path

    # ------------------------------------------------------------------
    # Decision-boundary probe
    # ------------------------------------------------------------------

    def decision_boundary_grid(self, x_range: Tuple[float, float] = (-0.5, 1.5),
                               y_range: Tuple[float, float] = (-0.5, 1.5),
                               n_grid: int = 80) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Evaluate the network output over a 2D input grid.

        Returns (X_mesh, Y_mesh, Z_pred) where Z_pred has shape (n, n).
        """
        xs = np.linspace(*x_range, n_grid)
        ys = np.linspace(*y_range, n_grid)
        X_mesh, Y_mesh = np.meshgrid(xs, ys, indexing="xy")
        flat = np.stack([X_mesh.ravel(), Y_mesh.ravel()], axis=1)
        out = self.predict(flat)
        if out.ndim > 1 and out.shape[1] == 1:
            out = out[:, 0]
        return X_mesh, Y_mesh, out.reshape(n_grid, n_grid)

    # ------------------------------------------------------------------
    # Generalisation / overfitting toy
    # ------------------------------------------------------------------

    @staticmethod
    def generate_noisy_dataset(n: int = 80, noise: float = 0.15, seed: int = 0
                               ) -> Tuple[np.ndarray, np.ndarray]:
        """Two interleaved Gaussian blobs in 2D for binary classification."""
        rng = np.random.default_rng(seed)
        n0, n1 = n // 2, n - n // 2
        c0 = rng.normal(loc=(0.0, 0.0), scale=noise, size=(n0, 2))
        c1 = rng.normal(loc=(1.0, 1.0), scale=noise, size=(n1, 2))
        X = np.vstack([c0, c1])
        y = np.concatenate([np.zeros(n0), np.ones(n1)])
        # Shuffle
        idx = rng.permutation(n)
        return X[idx], y[idx]

    @staticmethod
    def train_test_split(X: np.ndarray, y: np.ndarray,
                         test_fraction: float = 0.3, seed: int = 0
                         ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        rng = np.random.default_rng(seed)
        n = X.shape[0]
        idx = rng.permutation(n)
        n_test = max(1, int(round(test_fraction * n)))
        test_idx, train_idx = idx[:n_test], idx[n_test:]
        return X[train_idx], y[train_idx], X[test_idx], y[test_idx]

    def generalization_gap(self, X_train: np.ndarray, y_train: np.ndarray,
                           X_test: np.ndarray, y_test: np.ndarray) -> dict:
        """Loss on train vs test splits (positive gap => overfitting)."""
        L_tr = self._loss_value(self.predict(X_train), y_train)
        L_te = self._loss_value(self.predict(X_test),  y_test)
        return {
            "train_loss": float(L_tr),
            "test_loss":  float(L_te),
            "gap":        float(L_te - L_tr),
        }

    # ------------------------------------------------------------------
    # Substrate metadata
    # ------------------------------------------------------------------

    def substrate_signature(self) -> dict:
        return {
            "interpretation": (
                "training = strain-pattern relaxation under back-propagated "
                "error currents; gradient descent ~ viscous gradient flow; "
                "universality = completeness of substrate strain modes"
            ),
            "lambda_qcd_MeV": LAMBDA_QCD_MEV,
            "learning_rate": self.learning_rate,
            "loss": self.loss,
            "max_epochs": self.max_epochs,
            "tolerance": self.tolerance,
            "n_parameters": self.geometry.n_parameters(),
        }


__all__ = [
    "LAMBDA_QCD_MEV",
    "DEFAULT_LEARNING_RATE",
    "DEFAULT_MAX_EPOCHS",
    "DEFAULT_TOLERANCE",
    "DEFAULT_HIDDEN_WIDTH",
    "sigmoid",
    "sigmoid_prime",
    "tanh",
    "tanh_prime",
    "relu",
    "relu_prime",
    "boolean_dataset",
    "NeuralNetworkGeometry",
    "NeuralNetworkSimulator",
]
