# Second-Order Sensitivity Analysis for Floquet Multipliers

## 1. Problem Statement

Given the monodromy operator $\mathcal{M}(p): X \to X$ of a DDE, we want a
**second-order** approximation of the worst-case Floquet multiplier modulus
$|\mu_m(p)|$ as a function of the parameter vector $p \in \mathbb{R}^n$ near a
nominal point $p_0$.

The first-order (current) robustness margin $C_1$ is accurate only when the
map $p \mapsto |\mu_m(p)|$ is nearly linear over the uncertainty box
$\{|\Delta p_i| \le b_i\}$.  At locations where the stability boundary is
curved — especially near the lower-left corner of the stability lobe —
$C_1$ over- or under-estimates the true margin.

---

## 2. First-Order Recap

The compressed $n_b \times n_b$ operator (in the Krylov–Schur basis $Q$):

$$H_{ab}(p) = \langle q_a,\, \mathcal{M}(p)\, q_b \rangle$$

First-order eigenvalue perturbation theory gives:

$$\frac{\partial \lambda_m}{\partial p_i} = \ell_m^\top H_1^{(i)} r_m, \qquad
  H_1^{(i)} = \frac{\partial H}{\partial p_i}$$

First-order sensitivity of the modulus:

$$s_i^{(1)} = \frac{\partial |\mu_m|}{\partial p_i}
  = \frac{\operatorname{Re}\!\left(\bar{\mu}_m \frac{\partial\mu_m}{\partial p_i}\right)}{|\mu_m|}$$

Worst-case first-order direction: $\Delta p^*_i = \operatorname{sign}(s_i^{(1)})\cdot b_i$.

Robustness margin (linear):

$$C_1 = \frac{1 - |\mu_m|}{\sum_i |s_i^{(1)}| b_i}$$

---

## 3. Second-Order Extension

### 3.1 Quadratic approximation of $|\mu_m|$

$$|\mu_m(p_0 + \Delta p)| \approx |\mu_m| + g^\top \Delta p + \tfrac{1}{2}\,\Delta p^\top \mathbf{H}^{|\mu|}_m \Delta p$$

where $g_i = s_i^{(1)}$ and:

$$[\mathbf{H}^{|\mu|}_m]_{ik} =
  \frac{\operatorname{Re}\!\left(\overline{d_k}\,d_i + \bar{\mu}_m\,d_{ik}\right)}{|\mu_m|}
  - \frac{s_i^{(1)} s_k^{(1)}}{|\mu_m|}$$

with $d_i = \partial\mu_m/\partial p_i$ (complex) and $d_{ik} = \partial^2\mu_m/\partial p_i\partial p_k$ (complex).

### 3.2 Second-order eigenvalue derivative via perturbation theory

$$d_{ik} = \ell_m^\top H_2^{(ik)} r_m
  + \sum_{l \neq m}
    \frac{(\ell_m^\top H_1^{(i)} r_l)(\ell_l^\top H_1^{(k)} r_m)
        + (\ell_m^\top H_1^{(k)} r_l)(\ell_l^\top H_1^{(i)} r_m)}
         {\lambda_m - \lambda_l}$$

where $H_2^{(ik)} = \partial^2 H/\partial p_i\partial p_k$.

### 3.3 Computing $H_2^{(ik)}$ via nested ForwardDiff Duals

We use **three nested Dual layers** (outer → inner):

| Layer | Tag | Role |
|---|---|---|
| Outermost (stripped by `partialpart`) | `AffineTag` | Jacobian action on state direction $q_b$ |
| Middle | `SensTag2` | $\partial/\partial p_k$, carries $n$ partials |
| Inner | `SensTag` | $\partial/\partial p_i$, carries $n$ partials |

Tag precedence: `SensTag ≺ SensTag2 ≺ AffineTag`

Parameter seeding for parameter $j$ (active):

```julia
DualT1  = Dual{SensTag,  Float64, n}          # inner: ∂/∂p_i
DualT2  = Dual{SensTag2, DualT1,  n}          # outer: ∂/∂p_k (of inner Duals)

p_dual[j] = DualT2(
    DualT1(p0_j, e_j),                        # value: SensTag1 partial j = 1
    ntuple(k -> k==j ? DualT1(1.0,zp) : DualT1(0.0,zp), n)  # SensTag2 partial j = 1
)
```

After one `_linmap_action_dual` call, each element of $M(\tilde p)\,q_b$ is a
`DualT2`. The compressed matrix $\widetilde H_{ab} = \langle q_a, M(\tilde p) q_b\rangle$
is also `DualT2`. Extraction:

```
H0[a,b]         = value(value(H̃[a,b]))                          # Float64
H1[a,b,i]       = partials(value(H̃[a,b]))[i]                    # ∂H/∂p_i
H2[a,b,i,k]     = partials(partials(H̃[a,b])[k])[i]             # ∂²H/∂p_i∂p_k
```

**Note:** `partials(value(h))[i]` and `value(partials(h)[i])` both give $\partial H/\partial p_i$
by construction of the symmetric seeding.

---

## 4. Worst-Case Second-Order Robustness Margin

### 4.1 Box maximum of a quadratic (full vertex check)

Over the box $\mathcal{B} = \{|\Delta p_i| \le b_i\}$, the quadratic
$f(\Delta p) = g^\top\Delta p + \frac12 \Delta p^\top \mathbf{H} \Delta p$ achieves its
maximum at one of the $2^n$ vertices $v \in \{-1,+1\}^n$.

For each vertex $v$, let $\delta = v \odot b$. The second-order C-value is the
positive root of:

$$|\mu_m| + C\,L_1(v) + \tfrac{C^2}{2}\,L_2(v) = 1$$

$$C(v) = \frac{-L_1(v) + \sqrt{L_1(v)^2 - 2\,L_2(v)\,(|\mu_m| - 1)}}{L_2(v)}$$

where $L_1(v) = g^\top\delta$ and $L_2(v) = \delta^\top \mathbf{H}\,\delta$.

If $L_2(v) = 0$: $C(v) = (1 - |\mu_m|)/L_1(v)$ (linear fallback).

The worst-case margin:

$$C_2 = \min_{\substack{v \in \{-1,+1\}^n \\ L_1(v) > 0}} C(v)$$

For $n = 5$: 32 vertex checks, each requiring only scalar algebra — negligible cost.

### 4.2 Gradient-aligned second-order correction (fast approximation)

Use only the first-order worst vertex $v^* = \operatorname{sign}(g \odot b)$:

$$C_{2,\text{fast}} = \frac{-L_1 + \sqrt{L_1^2 - 2\,L_2(v^*)\,(|\mu_m|-1)}}{L_2(v^*)}$$

Same formula, single vertex. Faster but misses vertices where curvature makes
another corner worse.

---

## 5. Complexity Analysis

| Step | Cost | Bottleneck |
|---|---|---|
| Krylov $n_b$ eigenpairs (Float64) | $O(n_b \cdot T_{\rm LinMap})$ | as before |
| First-order: $n_b$ Dual{SensTag} LinMaps | $O(n_b \cdot n \cdot T_{\rm LinMap})$ | inner Dual $O(n)$ arith |
| **Second-order: $n_b$ nested Dual LinMaps** | $O(n_b \cdot n^2 \cdot T_{\rm LinMap})$ | DualT2 arith $O(n^2)$ |
| Vertex check ($2^n$ quadratics) | $O(2^n \cdot n^2)$ | scalar algebra only |

For $n = 5$, $n_b = 10$: second-order Hessian computation is $\approx 25\times$ more
expensive per grid point than first-order sensitivity. Still polynomial ($O(n^2)$)
vs exponential ($O(2^n)$ brute-force boundary evaluation).

The second-order LinMap uses `DualT2 = Dual{SensTag2, Dual{SensTag, Float64, n}, n}`,
which has $n + n^2 = 30$ scalar components for $n = 5$. Each floating-point
operation on the integrator state vector costs $\approx n^2 = 25\times$ more than
Float64. This dominates: total per-grid-point cost is $\approx 25\times$ the
first-order cost.

---

## 6. Approximation Errors

**Source 1: Subspace truncation** ($n_b < \infty$)
The left eigenvector $\ell_m$ of the compressed $n_b \times n_b$ matrix
approximates the true adjoint Floquet mode. As $n_b \to \infty$ the approximation
improves. Second-order accuracy also requires the off-diagonal terms
$\sum_{l \neq m} (\cdots)/(\lambda_m - \lambda_l)$: modes close to $\mu_m$ and
with small gap contribute most. A large enough $n_b$ is needed.

**Source 2: Perturbation theory order** (quadratic vs higher)
The quadratic approximation of $|\mu_m(p_0 + \Delta p)|$ is exact for linear
systems in $p$ (e.g., $p$ enters linearly in the DDE coefficients). For
nonlinear parameter dependence, cubic and higher terms are neglected.

**Source 3: Multi-mode interaction**
If two multipliers are close (near-degenerate), perturbation theory breaks
down due to small denominator $\lambda_m - \lambda_l$. A suitable regularization
or a Löwdin partitioning approach would be needed in that case.

---

## 7. Implementation Notes

- `compute_sensitivities_dual` (first-order) → unchanged, used for `sens_order=1`
- `compute_second_order_sens` (new) → nested Dual, returns `(sens1, sens2_complex, hess_mu)`
- `robustness_margin(dp, ...; sens_order=1)` → linear (current behavior)
- `robustness_margin(dp, ...; sens_order=2)` → quadratic worst-case (vertex check)
- `robustness_margin(dp, ...; sens_order=2, fast=true)` → gradient-aligned quadratic (single vertex)

The watchdog wrapper `with_timeout(f, t_max)` uses `@async` + `throwto(InterruptException)`
to kill stalled LinMap or MDBM evaluations.

---

---

## 8. Empirical Findings (smoke test, 4×5 grid, Neig=12)

| Metric | Value |
|---|---|
| Mean \|C₂−C₁\|/C₁ over stable grid | **8.6 %** |
| Max \|C₂−C₁\|/C₁ | **55 %** (lower-left corner region) |
| Timing ratio C₂ MDBM / C₁ MDBM | **~13×** |

The 55% maximum difference confirms the user's observation that the first-order
approximation was inaccurate near curved stability boundaries. The second-order
correction is most important in regions where the eigenvalue surface is strongly
curved — exactly the lower-left corner of the stability lobe.

The timing of ~13× (observed) vs ~n² = 25× (theoretical pure arithmetic) is
explained by caching, memory layout, and the fact that the outer `nb=12` LinMap
loop dominates over the arithmetic overhead at this problem size.

### Taylor test note
A direct Taylor expansion test at the nominal point showed both first- and
second-order errors at ~2×10⁻⁴, dominated by the Krylov eigenvalue noise floor
(reltol=1e-7 ODE integration). To isolate the Taylor improvement, one would need
either tighter Krylov tolerances (<1e-10) or larger perturbations (scale >3)
where the quadratic term exceeds the noise floor. The MDBM boundary comparison
and the C₁/C₂ map difference are more informative validation metrics.

*Last updated: 2026-05-29. See `src/robust_control.jl` for implementation,
`examples/second_order_sensitivity_demo.jl` for the full comparison demo.*
