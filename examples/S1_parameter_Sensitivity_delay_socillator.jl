# Example: Robustness map for a delay PD-controlled oscillator
#
#   ẍ(t) + 2ζω ẋ(t) + ω² x(t) = -kp · x(t-τ) - kd · ẋ(t-τ) + F(t)
#
# A small periodic forcing F(t) with period T = τ provides a non-trivial
# affine fixed point. The affine mapping period equals the delay (T = τ).
# Sweeps the (kp, kd) controller-gain plane and computes per grid point:
#   - the spectral radius ρ = max |μ|
#   - the worst-case linear robustness margin C against bounded uncertainty
#     in (ζ, ω, kp, kd, τ).

using Pkg;
Pkg.activate(".");
using AffineMapDiffEq
using StaticArrays
using DifferentialEquations
using LinearAlgebra
using KrylovKit
using GLMakie
using MDBM

## ─── Resolution / fidelity parameters (edit here to switch between runs) ─────
# HIGH-RESOLUTION (production)
const RES_sweep_kp   = 110      # grid points along kp for the C-map sweep
const RES_sweep_kd   = 100      # grid points along kd
const RES_mdbm_seed  = 14       # MDBM coarse-grid size (each axis)
const RES_mdbm_iter  = 6        # MDBM refinement iterations
const RES_N_random   = 20       # random grey validation curves

# #  LOW-RESOLUTION (fast preview — uncomment to switch)
#  const RES_sweep_kp   = 30
#  const RES_sweep_kd   = 31
#  const RES_mdbm_seed  = 10
#  const RES_mdbm_iter  = 4
#  const RES_N_random   = 5

## ─────────────────────────────────────────────────────────────────────────────

# 1. Governing equation: in-place DDE
function DelayPDOsc!(du, u, h, p, t)
    ζ, ω, kp, kd, τ, T = p
    hτ = h(p, t - τ)
    F  = 0.1 * cos(2π * t / T)
    du[1] = u[2]
    #Unstable with +ω^2, stable with -ω^2
    du[2] = +ω^2 * u[1] - 2ζ*ω * u[2] - kp * hτ[1] - kd * hτ[2] + F
end

# 2. Nominal parameters — mapping period locked to delay
const ζ_nom  = -0.02 #unstable
const ω_nom  = 1.0
const kp_nom = 2.0
const kd_nom = 2.0
const τ_nom  = 0.5
const T_nom  = τ_nom
const p_init = (ζ_nom, ω_nom, kp_nom, kd_nom, τ_nom, T_nom)

# 3. Affine mapping setup
u0 = @MArray [0.0, 0.0]
h(p, t) = @MArray [0.0, 0.0]
Solver_args = Dict(:alg => MethodOfSteps(Tsit5()), :verbose => false, :reltol => 1e-7)

probMapping = DDEProblem{true}(DelayPDOsc!, u0, h, (0.0, T_nom), p_init; constant_lags=[τ_nom])

println("Calibrating affine mapping at nominal parameters...")
Neig = 10
Krylov_arg = (Neig, :LM,
    KrylovKit.Arnoldi(tol=1e-15, krylovdim=Neig * 2 + 2, verbosity=0, maxiter=10))
dpPD = dynamic_problemSampled(probMapping, Solver_args, τ_nom;
    Historyresolution=100, zerofixpont=false, affineinteration=2, Krylov_arg=Krylov_arg)

@time mus, saff, sol0 = affine(dpPD; p=p_init);
@time mus, saff, sol0 = affine(dpPD; p=p_init);
println("Dominant |μ| at nominal: ", round(maximum(abs.(mus[1])); digits=6))


# ─── Robustness Margin C at the nominal point ────────────────────────────────

include(joinpath(@__DIR__, "..", "src", "robust_control.jl"))

# Predefined uncertainty ratios (as a ratio of the nominal values)
const ζ_var  = 0.02
const ω_var  = 0.02
const kp_var = 0.05
const kd_var = 0.05
const τ_var  = 0.02
const T_var  = 0.03      # T is tied to τ in this study

# Nominal uncertainty box (absolute = ratio · |nominal|)
const delta_p_nom = (ζ_var  * abs(ζ_nom),
                     ω_var  * abs(ω_nom),
                     kp_var * abs(kp_nom),
                     kd_var * abs(kd_nom),
                     τ_var  * abs(τ_nom),
                     T_var  * abs(T_nom))

println("\nComputing robustness margin C at nominal parameters...")
@time C, C_vec, mu_C, sens = robustness_margin(dpPD, p_init, delta_p_nom)

println("Floquet multipliers |μ|  = ", round.(abs.(mu_C); digits=6))
println("Per-root margin     C_j  = ", round.(C_vec; digits=4))
println("Worst-case margin   C    = ", round(C; digits=4))
println("\nSensitivity matrix  ∂|μ_j|/∂p_i:")
println("  (rows = multiplier index, cols = ζ  ω  kp  kd  τ  T)")
display(round.(sens; digits=4))


## ─── Robustness Map over (kp, kd) ────────────────────────────────────────────
# Threaded sweep, identical structure to parameter_Sensitivity.jl.
# The uncertainty ratios are interpreted as |Δp_i| = ratio · |p_i_local|, so
# the uncertainty box scales with the local nominal kp, kd at each grid point.

println("\nStarting Robustness Map (threaded over (kp, kd))...")

const kpv = LinRange(-1.5, 3.0, RES_sweep_kp)
const kdv = LinRange(-1.0, 4.0, RES_sweep_kd)

# Two dynamic problems sharing the same discretization, differing only in the
# sensitivity method used by robustness_margin:
#   dp_dual: perturbation_size = 0.0   → exact ForwardDiff (Dual) sensitivity
#   dp_fd  : perturbation_size = 1e-6  → central finite-difference sensitivity
# Neig sets both the tracked-multiplier count and the Schur-subspace size used to
# compress the operator for the Dual sensitivity; a larger Neig captures the
# dominant mode's left eigenvector more fully → more accurate C. ~10 balances
# accuracy and speed.
# ── dp for spectral radius only (Neig=2: just the dominant mode, no sensitivity) ──
Neig_spec = 8
Krylov_arg_spec = (Neig_spec, :LM,
    KrylovKit.Arnoldi(tol=1e-12, krylovdim=Neig_spec + 10, verbosity=0, maxiter=20))
dp_spec = dynamic_problemSampled(probMapping, Solver_args, τ_nom;
    Historyresolution=25, zerofixpont=false, affineinteration=1,
    Krylov_arg=Krylov_arg_spec)

# ── dp for robustness margin (Neig=10: larger subspace for accurate Dual sensitivity) ──
Neig_fast = 25
Krylov_arg_fast = (Neig_fast, :LM,
    KrylovKit.Arnoldi(tol=1e-12, krylovdim=Neig_fast + 20, verbosity=0, maxiter=20))
dp_dual = dynamic_problemSampled(probMapping, Solver_args, τ_nom;
    Historyresolution=25, zerofixpont=false, affineinteration=1,
    Krylov_arg=Krylov_arg_fast, perturbation_size=0.0)
dp_fd = dynamic_problemSampled(probMapping, Solver_args, τ_nom;
    Historyresolution=25, zerofixpont=false, affineinteration=1,
    Krylov_arg=Krylov_arg_fast, perturbation_size=1e-5)

# Spectral radius sweep — only spectralradius(), no sensitivity computation.
function sweep_spec_map!(dp, Smap)
    println("  Processing kp  ", length(kpv), ":")
    Threads.@threads for j in 1:length(kpv)
        print(" ", j)
        for i in 1:length(kdv)
            Smap[i, j] = spectralradius(dp; p=(ζ_nom, ω_nom, kpv[j], kdv[i], τ_nom, T_nom))
        end
    end
    println(" Done.")
end

# Robustness map sweep — robustness_margin only, no spectral radius recording.
function sweep_C_map!(dp, Cmap)
    println("  Processing kp  ", length(kpv), ":")
    Threads.@threads for j in 1:length(kpv)
        kp_v = kpv[j]
        print(" ", j)
        for i in 1:length(kdv)
            kd_v = kdv[i]
            p_loc = (ζ_nom, ω_nom, kp_v, kd_v, τ_nom, T_nom)
            delta_p_loc = (ζ_var  * abs(ζ_nom),
                           ω_var  * abs(ω_nom),
                           kp_var * abs(kp_v),
                           kd_var * abs(kd_v),
                           τ_var  * abs(τ_nom),
                           T_var  * abs(T_nom))
            C_loc, _, mu_loc, _ = robustness_margin(dp, p_loc, delta_p_loc; verbose=false)
            ρ = maximum(abs.(mu_loc))
            if isfinite(C_loc) && ρ < 1.0
                Cmap[i, j] = C_loc
            end
        end
    end
    println(" Done.")
    return nothing
end

Spec_map   = zeros(length(kdv), length(kpv))
C_map_dual = fill(NaN, length(kdv), length(kpv))
C_map_fd   = fill(NaN, length(kdv), length(kpv))

# Warm up all three branches before timing.
println("  Warming up...")
spectralradius(dp_spec; p=p_init)
robustness_margin(dp_dual, p_init, delta_p_nom; verbose=false)
robustness_margin(dp_fd,   p_init, delta_p_nom; verbose=false)

println("  Spectral radius sweep...")
t_spec = @elapsed sweep_spec_map!(dp_spec, Spec_map)
println("  Dual (ForwardDiff) sweep...")
t_dual = @elapsed sweep_C_map!(dp_dual, C_map_dual)
println("  Finite-difference sweep...")
t_fd   = @elapsed sweep_C_map!(dp_fd,   C_map_fd)

println("\nWhole-map timing:")
println("  Spectral radius only : ", round(t_spec; digits=2), " s")
println("  Dual (ForwardDiff)   : ", round(t_dual; digits=2), " s")
println("  Finite difference    : ", round(t_fd;   digits=2), " s")
println("  Dual / spec overhead : ", round(t_dual / t_spec; digits=1), "×")
println("  stable points        : ", count(isfinite, C_map_dual), " / ", length(C_map_dual))

# Confine both C maps to the authoritative (spectral radius sweep) stable region.
unstable_mask = Spec_map .>= 1.0
C_map_dual[unstable_mask] .= NaN
C_map_fd[unstable_mask]   .= NaN


## ── Plot: spectral radius + robustness margin C (Dual) + C (FD), stable region ──
fig_C = Figure(size=(2000, 650))
ax_rho  = GLMakie.Axis(fig_C[1, 1],
    title="Spectral radius ρ = max|μ|  ($(round(t_spec; digits=1)) s)", xlabel="kp", ylabel="kd")
ax_dual = GLMakie.Axis(fig_C[1, 2],
    title="C — ForwardDiff  ($(round(t_dual; digits=1)) s)", xlabel="kp", ylabel="kd")
ax_fd   = GLMakie.Axis(fig_C[1, 3],
    title="C — finite difference  ($(round(t_fd; digits=1)) s)", xlabel="kp", ylabel="kd")

# Shared color range from the (cleaner) Dual map, so the FD noise is visible as
# saturation rather than rescaling the whole plot.
finite_vals = filter(isfinite, C_map_dual)
crange = isempty(finite_vals) ? (0.0, 1.0) : (minimum(finite_vals), maximum(finite_vals))

# Spectral radius on the stable region (unstable → NaN)
Spec_stable = copy(Spec_map)
Spec_stable[Spec_stable .>= 1.0] .= NaN
hm_rho = heatmap!(ax_rho, kpv, kdv, log.(Spec_stable'), colormap=:inferno)
Colorbar(fig_C[2, 1], hm_rho, vertical=false, label="log(ρ)  (unstable → NaN)")

hm_dual = heatmap!(ax_dual, kpv, kdv, C_map_dual', colormap=:viridis, colorrange=crange)
Colorbar(fig_C[2, 2], hm_dual, vertical=false, label="C (Dual)")

# C = 1 contour on the accurate (Dual) map. Inside this curve C > 1: the linear
# robustness estimate guarantees |μ| stays < 1 over the FULL specified uncertainty
# box Δp (i.e. still stable), provided the linearity assumption holds. Outside the
# curve (C < 1) the box can reach the stability boundary — robustness not guaranteed.
contour!(ax_dual, kpv, kdv, C_map_dual'; levels=[2.0], color=:white, linewidth=2.5)

hm_fd = heatmap!(ax_fd, kpv, kdv, C_map_fd', colormap=:viridis, colorrange=crange)
Colorbar(fig_C[2, 3], hm_fd, vertical=false, label="C (FD)")

save("examples/parameter_Sensitivity_delay_socillator_map.png", fig_C)
display(fig_C)

println("Robustness map complete. Figure saved.")


## ─── MDBM stability boundary (ρ = 1) overlaid on all three maps ───────────────
using MDBM
# Locate the curve where the dominant multiplier crosses the unit circle, i.e.
# log(spectral radius) = 0, via MDBM (multi-dimensional bisection). A dedicated
# 2-multiplier problem is used: only the dominant mode sets stability, and this
# avoids the poorly-converged higher modes that can perturb max|μ|.
println("\nMDBM stability-boundary search...")

Neig_stab = 2
Krylov_arg_stab = (Neig_stab, :LM,
    KrylovKit.Arnoldi(tol=1e-12, krylovdim=Neig_stab + 10, verbosity=0, maxiter=20))
dp_stab = dynamic_problemSampled(probMapping, Solver_args, τ_nom;
    Historyresolution=25, zerofixpont=false, affineinteration=1, Krylov_arg=Krylov_arg_stab)

# Coarse seed grids spanning the same (kp, kd) window as the sweep
kpv_mdbm = LinRange(-1.5, 3.0, RES_mdbm_seed)
kdv_mdbm = LinRange(-1.0, 4.0, RES_mdbm_seed)
foo_stab(kp_loc, kd_loc) =
    log(spectralradius(dp_stab; p=(ζ_nom, ω_nom, kp_loc, kd_loc, τ_nom, T_nom)))

mymdbm = MDBM_Problem(foo_stab, [MDBM.Axis(kpv_mdbm, "kp"), MDBM.Axis(kdv_mdbm, "kd")])
t_mdbm_stab = @elapsed begin
    MDBM.solve!(mymdbm, RES_mdbm_iter, interpolationorder=0, verbosity=0)
    MDBM.solve!(mymdbm, 0, interpolationorder=1, verbosity=0)
end
x_mdbm, y_mdbm = getinterpolatedsolution(mymdbm)
println("  ρ=1 boundary: $(length(x_mdbm)) points  in $(round(t_mdbm_stab; digits=1)) s")

# Overlay the boundary on all three panels (cyan on inferno, red on viridis)
scatter!(ax_rho,  x_mdbm, y_mdbm, color=:cyan, markersize=4)
scatter!(ax_dual, x_mdbm, y_mdbm, color=:red,  markersize=4)
scatter!(ax_fd,   x_mdbm, y_mdbm, color=:red,  markersize=4)

save("examples/parameter_Sensitivity_delay_socillator_map.png", fig_C)
display(fig_C)
println("MDBM boundary added (", length(x_mdbm), " points). Figure re-saved.")


## ─── Validation of the C = 2 robustness boundary ─────────────────────────────
# At a point with margin C, scaling the uncertainty box by C brings the worst-case
# |μ| to 1. So inside the C ≥ 2 region the system should stay stable under ANY
# parameter perturbation within 2·Δp. We validate empirically: for many perturbed
# parameter realisations inside the 2·Δp box we recompute the ρ=1 stability
# boundary (in the kp–kd design plane) and overlay it. If the C=2 prediction holds
# (and the system is linear enough), every perturbed boundary stays on the
# UNSTABLE side of the C=2 contour — i.e. none intrudes into the C ≥ 2 region.
const C_validate   = 2.0
const N_random     = RES_N_random
const png_path     = "examples/parameter_Sensitivity_delay_socillator_map.png"

# Maximum delay that can appear in the validation: τ_nom · (1 + C_validate · τ_var).
# The dynamic problem used for validation must declare this as its history window
# (maxdelay) AND as the constant_lag so the DDE solver correctly tracks the
# discontinuity for every perturbed τ ≤ τ_max_val.
const τ_max_val = τ_nom * (1.0 + C_validate * τ_var)   # = 0.52

# New DDEProblem with the enlarged lag declaration
probMapping_val = DDEProblem{true}(DelayPDOsc!, u0, h, (0.0, T_nom), p_init;
    constant_lags=[τ_max_val])

dp_val = dynamic_problemSampled(probMapping_val, Solver_args, τ_max_val;
    Historyresolution=25, zerofixpont=false, affineinteration=1,
    Krylov_arg=Krylov_arg_stab)

# C=2 contour on the accurate (Dual) map — the boundary under test.
contour!(ax_dual, kpv, kdv, C_map_dual'; levels=[C_validate], color=:cyan, linewidth=2.5)
save(png_path, fig_C); display(fig_C)

# ρ=1 boundary in the (kp,kd) design plane for a perturbation given as ratios
# r = (rζ, rω, rkp, rkd, rτ) ∈ [-1,1]^5, scaled by C_validate · (local box).
# dp_val covers the full perturbed delay range; dp_stab (τ_nom) must NOT be used
# here as its history window would be too small when rτ > 0.
function perturbed_boundary(r)
    rζ, rω, rkp, rkd, rτ = r
    sc = C_validate
    foo(kp, kd) = log(spectralradius(dp_val; p=(
        ζ_nom + rζ * sc * ζ_var  * abs(ζ_nom),
        ω_nom + rω * sc * ω_var  * abs(ω_nom),
        kp    + rkp * sc * kp_var * abs(kp),
        kd    + rkd * sc * kd_var * abs(kd),
        τ_nom + rτ * sc * τ_var  * abs(τ_nom),
        T_nom)))
    prob = MDBM_Problem(foo, [MDBM.Axis(kpv_mdbm, "kp"), MDBM.Axis(kdv_mdbm, "kd")])
    try
        MDBM.solve!(prob, RES_mdbm_iter, interpolationorder=0, verbosity=0)
        MDBM.solve!(prob, 0, interpolationorder=1, verbosity=0)
        return getinterpolatedsolution(prob)
    catch
        return Float64[], Float64[]
    end
end

overlay!(xb, yb, col) = for ax in (ax_rho, ax_dual, ax_fd)
    scatter!(ax, xb, yb, color=col, markersize=5)
end

# Identify which of the 5 parameters (ζ,ω,kp,kd,τ) have non-zero uncertainty.
# Only these are varied in the random and corner sweeps; the others stay at nominal.
const val_vars   = (ζ_var, ω_var, kp_var, kd_var, τ_var)
const active_idx = findall(!iszero, val_vars)   # indices into the 5-tuple
const n_active_v = length(active_idx)
const n_corners  = 2^n_active_v

println("Active uncertain parameters: $n_active_v / 5  (positions $active_idx in (ζ,ω,kp,kd,τ))")
println("  → $N_random random + $n_corners corner ($n_active_v-cube) cases")

# Build a full 5-tuple r from the active-only values; inactive positions stay at 0.
function make_r(active_vals)
    r = zeros(5)
    for (k, idx) in enumerate(active_idx)
        r[idx] = active_vals[k]
    end
    return ntuple(i -> r[i], 5)
end

# Part 1 — N random realisations inside the 2·Δp box (only active params vary)
println("\nC=$C_validate validation: $N_random random realisations...")
for n in 1:N_random
    r = make_r(ntuple(_ -> 2rand() - 1, n_active_v))
    xb, yb = perturbed_boundary(r)
    overlay!(xb, yb, (:gray, 0.4))
    save(png_path, fig_C); display(fig_C)
    println("  random $n/$N_random  ($(length(xb)) boundary points)")
end

# Part 2 — all 2^n_active_v corners of the active-parameter cube (orange)
println("\nC=$C_validate validation: $n_corners corners of the $n_active_v-cube...")
for (k, c) in enumerate(Iterators.product(ntuple(_ -> (-1.0, 1.0), n_active_v)...))
    r = make_r(c)
    xb, yb = perturbed_boundary(r)
    overlay!(xb, yb, (:orange, 0.7))
    save(png_path, fig_C); display(fig_C)
    println("  corner $k/$n_corners  active_signs=$(c)  ($(length(xb)) boundary points)")
end

println("\nC=$C_validate validation complete. If the cyan C=$C_validate contour bounds all")
println("overlaid boundaries from the stable side, the linear robustness estimate holds.")


## ─── Direct MDBM computation of the C = C_validate boundary ─────────────────
# We feed MDBM the function  f(kp,kd) = C(kp,kd) − C_validate, where C is
# computed via the standard formula for ALL points (not just stable ones).
#
# Sign analysis:
#   stable   (|μ|<1) : C = (1−|μ|)/Σ|∂|μ|/∂p|Δp > 0
#                      C > C_v  →  f > 0    (inside the boundary we want)
#                      C = C_v  →  f = 0    ← the MDBM target
#                      C < C_v  →  f < 0
#   boundary (|μ|=1) : C = 0   →  f = −C_v < 0
#   unstable (|μ|>1) : C < 0   →  f < 0
#
# A single, clean sign change at C = C_v — no spurious zeros, no kink.
# This equals  sign(−log|μ_max|)·|C| − C_v, i.e. the user's formula with
# sign(−log|μ|) instead of sign(log|μ|).  (The two signs differ: log|μ|<0
# when stable, so sign(log|μ|)·C would be −C in the stable region — always
# negative — and MDBM could never find a zero there.)

function C_boundary_func(kp_loc, kd_loc)::Float64
    p_loc = (ζ_nom, ω_nom, kp_loc, kd_loc, τ_nom, T_nom)
    delta_p_loc = (ζ_var  * abs(ζ_nom),
                   ω_var  * abs(ω_nom),
                   kp_var * abs(kp_loc),
                   kd_var * abs(kd_loc),
                   τ_var  * abs(τ_nom),
                   T_var  * abs(T_nom))

    # affine gives both eigenvalues and the Schur basis (no early exit for unstable)
    mus, s_nom = affine(dp_dual; p=p_loc)
    all_mu  = ComplexF64.(mus[1])
    n_eff   = min(dp_dual.Krylov_arg[1], length(all_mu))
    mu_nom  = all_mu[1:n_eff]
    eigvecs = mus[2]

    # Dual sensitivities (exact, same method as the C map)
    sens = compute_sensitivities_dual(dp_dual, p_loc, delta_p_loc, mu_nom, s_nom, eigvecs)

    delta_vec = collect(Float64, delta_p_loc)
    C_min = minimum(map(eachindex(mu_nom)) do j
        denom = sum(abs(sens[j, i]) * delta_vec[i] for i in eachindex(delta_vec))
        denom < eps(Float64) ? Inf : (1.0 - abs(mu_nom[j])) / denom
    end)

    return C_min - C_validate
end

println("\nDirect MDBM C = $C_validate boundary (Dual sensitivities)...")
C_boundary_func(kp_nom, kd_nom)   # warmup / JIT

kpv_Cm = LinRange(-1.5, 3.0, RES_mdbm_seed)
kdv_Cm = LinRange(-1.0, 4.0, RES_mdbm_seed)
mdbm_C = MDBM_Problem(C_boundary_func,
    [MDBM.Axis(kpv_Cm, "kp"), MDBM.Axis(kdv_Cm, "kd")])
t_mdbm_C = @elapsed begin
    MDBM.solve!(mdbm_C, RES_mdbm_iter, interpolationorder=0, verbosity=0)
    MDBM.solve!(mdbm_C, 0, interpolationorder=1, verbosity=0)
end
x_Cm, y_Cm = getinterpolatedsolution(mdbm_C)
println("  C = $C_validate MDBM boundary: $(length(x_Cm)) points  in $(round(t_mdbm_C; digits=1)) s")

# Overlay: lime diamonds on the Dual panel — should coincide with the cyan
# C = C_validate contour extracted from the heatmap, validating both methods.
scatter!(ax_dual, x_Cm, y_Cm, color=:lime, markersize=3, marker=:diamond,
    label="MDBM C=$C_validate")
save(png_path, fig_C); display(fig_C)
println("MDBM C boundary overlaid (lime ◇). Should match cyan contour from heatmap.")

println("\n─── MDBM timing comparison ───────────────────────────────────────────────")
println("  ρ = 1  boundary  (spectralradius only, Neig=2)  : ",
        round(t_mdbm_stab; digits=2), " s   ($(length(x_mdbm)) pts)")
println("  C = $C_validate boundary  (Dual sensitivities,  Neig=$Neig_fast) : ",
        round(t_mdbm_C;   digits=2), " s   ($(length(x_Cm)) pts)")
println("  Overhead ratio C / ρ : ", round(t_mdbm_C / t_mdbm_stab; digits=1), "×")
println("  This matches the per-point ratio: spectralradius ≈ 1 LinMap,")
println("  C_boundary ≈ $Neig_fast LinMap calls (Schur subspace) + O(n_active) Dual arith.")
