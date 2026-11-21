"""chem_factor: shared numpy / JAX pieces for the DF-THC factorization work (FeMoCo 76o, official verifier data).

The starter kit's `dfthc` trainer is not on this machine, so the objective is rebuilt from the verifier's own
reference loader: error_2 in the packed space of orbital pairs,

    A[(r, c), p] = sum_b w[r, b, c] u[r, b, i_p] u[r, b, j_p]          p = (i_p >= j_p), P = 2926
    error_2      = || sqrt(m m^T) * (A^T A - T) ||_F ,  T = tei + (S (x) I + I (x) S) / 2

with m = 1 on diagonal pairs and 2 elsewhere. It equals the verifier's error_2 (checked in `selfcheck`).
"""
from __future__ import annotations

import os
import sys

import numpy as np

DATA = os.path.expanduser("~/qcc/data")
os.environ.setdefault("DFTHC_DATA", DATA)
PIPE = os.path.expanduser("~/qcc/src/marculera_pipeline/dfthc_refit")
for _p in (os.path.join(PIPE, "experiments"), PIPE, os.path.join(DATA, "verifiers"), os.path.join(DATA, "kickmix")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

N = 76
UPPER_ENERGY = -22140.27335
GAP_FLOOR = 3.914898


def log_depth_tree(n: int = N) -> np.ndarray:
    edges, level, depth = [], [0], 0
    while len(edges) < n - 1:
        nxt = []
        for label in level:
            nxt.append(label)
            right = label + (1 << depth)
            if right < n:
                edges.append((label, right))
                nxt.append(right)
        level, depth = nxt, depth + 1
    return np.asarray(edges, dtype=np.int64)


def levels(edges: np.ndarray) -> list[np.ndarray]:
    """Runs of consecutive edges touching disjoint orbitals (one vectorized update each)."""
    runs, current, used = [], [], set()
    for t, (p, q) in enumerate(edges):
        if p in used or q in used:
            runs.append(np.array(current))
            current, used = [], set()
        current.append(t)
        used.update((int(p), int(q)))
    return runs + [np.array(current)]


def decode_np(theta_over_pi: np.ndarray, edges: np.ndarray) -> np.ndarray:
    th = np.pi * np.asarray(theta_over_pi, np.float64)
    v = np.zeros((th.shape[0], edges.shape[0] + 1))
    v[:, 0] = 1.0
    for t, (p, c) in enumerate(edges):
        a, b = v[:, p].copy(), v[:, c].copy()
        v[:, p] = np.cos(th[:, t]) * a - np.sin(th[:, t]) * b
        v[:, c] = np.sin(th[:, t]) * a + np.cos(th[:, t]) * b
    return v


def signed_angles(vectors: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """theta/pi in [-1/2, 1/2] reproducing each row up to a global sign (pipeline's encode.signed_angles)."""
    a = np.array(vectors, dtype=np.float64)
    a /= np.linalg.norm(a, axis=1, keepdims=True)
    theta = np.zeros((len(a), len(edges)))
    for t in range(len(edges) - 1, -1, -1):
        p, q = edges[t]
        s = np.where(a[:, p] > 0, 1.0, np.where(a[:, p] < 0, -1.0, np.where(a[:, q] > 0, -1.0, 1.0)))
        theta[:, t] = np.arctan2(s * a[:, q], s * a[:, p])
        a[:, p] = s * np.hypot(a[:, p], a[:, q])
        a[:, q] = 0.0
    return theta / np.pi


class Reference:
    """Reference integrals, the hosts' fixed shift and the packed two-body target."""

    def __init__(self):
        from dfthcsos.scoreinfra.reference import check_fcidump, load_reference
        from dfthcsos.scoreinfra.dfthc_math import symmetric_from_lower
        ref = load_reference(check_fcidump())
        self.ref = ref
        with np.load(os.path.join(DATA, "verifiers", "dfthcsos", "FeMoCo76o113e.fixed_shifts.npz")) as z:
            self.shifts = {k: np.array(z[k]) for k in ("H1_I", "Hsym_lower")}
        self.S = np.asarray(symmetric_from_lower(self.shifts["Hsym_lower"], N), np.float64)
        self.beta = float(self.shifts["H1_I"]) - 0.5 * float(np.trace(self.S))
        self.eta = int(ref["electrons"])
        tei = np.asarray(ref["tei"], np.float64)
        self.h1_thc = np.asarray(ref["oei"], np.float64) - 0.5 * np.einsum("prrq->pq", tei)
        self.h1_base = self.h1_thc - 0.5 * self.eta * self.S + self.beta * np.eye(N)
        self.energy = float(UPPER_ENERGY - ref["ecore"] + self.eta * self.beta)
        i, j = np.tril_indices(N)
        self.i, self.j = i, j
        T = tei[i, j][:, i, j].copy()                                  # (P, P)
        eye = (i == j).astype(np.float64)
        Sp = self.S[i, j]
        T += 0.5 * (np.outer(Sp, eye) + np.outer(eye, Sp))
        self.T = T
        m = np.where(i == j, 1.0, 2.0)
        self.M = np.outer(m, m)
        self.T2 = float(np.sum(self.M * T * T))

    def error2_np(self, u: np.ndarray, w: np.ndarray) -> float:
        R, B, C = w.shape
        u = u.reshape(R, B, N)
        u = u / np.linalg.norm(u, axis=2, keepdims=True)
        pairs = u[:, :, self.i] * u[:, :, self.j]
        A = np.einsum("rbc,rbp->rcp", w, pairs).reshape(R * C, -1)
        H = A.T @ A - self.T
        return float(np.sqrt(np.sum(self.M * H * H)))

    def official(self, words: np.ndarray, w: np.ndarray, W: np.ndarray, edges: np.ndarray, brot: int = 8) -> dict:
        """The verifier's metrics for a full word table (one-body rows first), weights w and offsets W."""
        from dfthcsos.scoreinfra.dfthc_math import evaluate_dfthc_angle_words
        m = evaluate_dfthc_angle_words(words, w, W, self.S, self.beta, self.ref, brot=brot, tree_edges=edges,
                                       energy_upper=UPPER_ENERGY)
        return {k: float(m[k]) for k in ("error_1", "error_2", "Lambda", "E_gap", "lambda_eff")}


def formula_costs(R: int, B: int, C: int, brot: int = 8, J: int = 5, angles_per_pass: float = 0.0, angles_temp: int = 77) -> dict:
    """T_total and Q_total of the verifier's cost model for a shape, given the lookup Toffolis per pass."""
    from dfthcsos.scoreinfra.costs_prepare import SpinFreeData_V1
    s = SpinFreeData_V1(R=R, B=B, C=C, N=N, chi=10, beta=brot)
    outer, inner, giv = s.outer_state_prep(), s.inner_state_prep(), s.givens_rotations()
    qo, qi = s.outer_state_prep_qubit_cost(), s.inner_state_prep_qubit_cost()
    Q = qo["persistent"] + qi["persistent"] + J * brot + max(qo["temporary"], qi["temporary"], angles_temp)
    return {"outer": outer, "inner": inner, "givens": giv, "heuristic": outer + inner + giv,
            "T": outer + inner + giv + 2 * angles_per_pass, "Q": Q, "rows": N + R * B, "blocks": -(-(N + R * B) // 16)}


# ----------------------------------------------------------------------------------------------- JAX side

def jax_setup(x64: bool = True):
    import jax
    jax.config.update("jax_enable_x64", bool(x64))
    jax.config.update("jax_default_matmul_precision", "highest")
    import jax.numpy as jnp
    return jax, jnp


def make_decode(edges: np.ndarray):
    import jax.numpy as jnp
    lv = levels(edges)

    def decode(theta):
        c, s = jnp.cos(jnp.pi * theta), jnp.sin(jnp.pi * theta)
        u = jnp.zeros((theta.shape[0], N), theta.dtype).at[:, 0].set(1.0)
        for ts in lv:
            p, q = edges[ts, 0], edges[ts, 1]
            up, uq = u[:, p], u[:, q]
            u = u.at[:, p].set(c[:, ts] * up - s[:, ts] * uq).at[:, q].set(s[:, ts] * up + c[:, ts] * uq)
        return u
    return decode


def make_error2(refc: Reference, dtype):
    """error2(u (R,B,N) unit rows, w (R,B,C)) and residual(u, w) (weighted, P x P) as JAX functions."""
    import jax.numpy as jnp
    T = jnp.asarray(refc.T, dtype)
    sM = jnp.asarray(np.sqrt(refc.M), dtype)
    i, j = jnp.asarray(refc.i), jnp.asarray(refc.j)

    def residual(u, w):
        pairs = u[:, :, i] * u[:, :, j]
        A = jnp.einsum("rbc,rbp->rcp", w, pairs).reshape(-1, T.shape[0])
        return sM * (A.T @ A - T)

    def error2(u, w):
        r = residual(u, w)
        return jnp.sqrt(jnp.sum(r * r))
    return error2, residual


def d_penalty(w):
    """Sign mixing D(w) = 1/2 sum_rc ((sum_b |w|)^2 - (sum_b w)^2): the part of Lambda the fit controls."""
    import jax.numpy as jnp
    return 0.5 * jnp.sum(jnp.sum(jnp.abs(w), axis=1) ** 2 - jnp.sum(w, axis=1) ** 2)


def selfcheck():
    """error2_np against the official verifier on the public table."""
    pub = os.path.expanduser("~/qcc/work/pub/dfthcsos/submission.npz")
    with np.load(pub) as z:
        words, w, h2i, edges = z["angle_words"], z["W"], z["H2_I"], z["tree_edges"]
    refc = Reference()
    u = decode_np(words[N:] / 128.0, edges)
    e2 = refc.error2_np(u, w)
    m = refc.official(words, w, h2i - w.sum(axis=1), edges)
    print("packed error_2", e2, "official", m)
    print(formula_costs(*w.shape, angles_per_pass=4133))


if __name__ == "__main__":
    selfcheck()
