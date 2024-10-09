5 + 5
using Interpolations
using StaticArrays
using DifferentialEquations
using Plots
plotly()
gr()
#using Contour
using PlotlyJS
using MAT
using LinearAlgebra
using LaTeXStrings
using FFTW
using ForwardDiff
using Revise
using KrylovKit
using Base.Threads

includet("src\\DDE_mapping_types.jl")
includet("src\\DDE_mapping_functions.jl")

Base.:+(a::SVector, b::Bool) = a .+ b

# ----------------------- Parameters -----------------------------

w = 0.03 / 1e3;

ϕshift = deg2rad(45.0);

Ωrel = 170 / 1.55;
RVA   = 0.1;
RVF   = 1 / 2;

m  = 2.9214;
k  = 5.8995e+6;
c  = 40.4073;
f0 = 0.0001;

ωn = sqrt(k / m);
ζ  = c / (2 * ωn * m);
K  = 431.25e+6;
Tn = 2pi / ωn;

T = maximum([2pi / (Ωrel * RVF), 2pi / Ωrel]);
p = [ζ, ωn, k, w, K, Ωrel, RVA, RVF, f0, ϕshift];

Ω(t)         = Ωrel*(1.0 + RVA * sin(t * Ωrel * RVF))
ϕ(t, ϕshift) = Ωrel * t + ϕshift - RVA / RVF * (cos(RVF * Ωrel * t) - 1);
g(ϕ)        = all([deg2rad(0) < mod(ϕ, 2pi), mod(ϕ, 2pi) < deg2rad(30)]) ? 0 : 1

function plot_phase()
    ϕ_vec = ϕ.(range(0,T,100), deg2rad(180))
    g_vec  = g.(ϕ_vec)
    Ω_vec  = Ω.(range(0,T,100))
    
    plot((ϕ_vec), g_vec, xlabel="ϕ", ylabel="g", label="g(ϕ)")
    plot!((ϕ_vec), Ω_vec./Ωrel, label="Ω(ϕ)/Ωrel")
end

@show lobenumber_range = ωn ./ (Ωrel .* (1.0 .+ [RVA, -RVA]))

# -------------------------------------------------------------------------------------------------------------
# ------------------------------------- Testing affine mapping ------------------------------------------------
# -------------------------------------------------------------------------------------------------------------

algMoS      = MethodOfSteps(BS3())
Krylov_arg  = (2, :LM, KrylovKit.Arnoldi(tol=1e-10, krylovdim=5, verbosity=0));
Solver_args = Dict(:alg => MethodOfSteps(BS3()), :verbose => false, :reltol => 1e-5, :maxiters => 1e7,)

u0       = SA[0.0, 0.0];
h(p, t)  = SA[0.0; 0.0]
τmax     = 2pi / Ωrel / (1.0 - RVA)
dtapprox = Tn / 30;

@show Nstep = Int(ceil(τmax / dtapprox))

function delayed_turning_STATIC_SSV(u, h, p, t)
    ζ, ωn, k, w, K, Ωrel, RVA, RVF, f0, ϕshift = p
    
    T0  = 2pi / Ωrel
    τ = T0 / (1.0 + RVA * sin(t * Ωrel * RVF))
    h0  = f0 * τ / T0

    dx  = u[2]
    ddx = ωn^2 * (g(ϕ(t, ϕshift)) * K * w / k * (h(p, t - τ)[1] - u[1] - h0) - u[1]) - 2 * ζ * ωn * u[2]

    SA[dx, ddx]
end
probTurningSSV = DDEProblem(delayed_turning_STATIC_SSV, u0, h, [0, T], p)
dpdp = dynamic_problemSampled(probTurningSSV, Solver_args, τmax, T;
    Historyresolution = Nstep, 
    zerofixpont       = false, 
    affineinteration  = 1,
    Krylov_arg        = Krylov_arg);

@time muaff, s0aff, solfix = affine(dpdp; p=[ζ, ωn, k, w, K, Ωrel, RVA, RVF, f0, ϕshift]);

fnan(y) = y == 0 ? NaN : y
plot(solfix.t, real.(getindex.(solfix.u, 1)))
plot!(solfix.t, fnan.(g.(ϕ.(solfix.t, ϕshift)) .* real.(getindex.(solfix.u, 1))))















# -------------------------------------------------------------------------------------------------------------
# ------------------------------------- p2p(w,ϕshift) surface and contour plot with heatmap -------------------
# -------------------------------------------------------------------------------------------------------------
RVF = 1
RVA = 0.1
Ωrel_vec   = 230

ϕshift_vec  = range(0.0, 4*pi,10)
w_vec       = range(0.01, 2.5, 10) / 1e3

Aaff     = zeros(size(ϕshift_vec, 1), size(w_vec, 1))
Spek_aff = zeros(size(ϕshift_vec, 1), size(w_vec, 1))
p2p      = zeros(size(ϕshift_vec, 1), size(w_vec, 1))
p2p_stab = zeros(size(ϕshift_vec, 1), size(w_vec, 1))


println("---------- Starting brute force ---------------")
@time for (ii,ϕshift) in enumerate(ϕshift_vec), (jj,w) in enumerate(w_vec)
    
    T     = maximum([2pi / (Ωrel * RVF), 2pi / Ωrel]);
    τmax  = 2pi / Ωrel / (1.0 - RVA) * 1.0
    Nstep = Int(ceil(τmax / dtapprox))

    function fooTurningSSV(Ωrel, w)

        dpdploc = dynamic_problemSampled(probTurningSSV, Solver_args, τmax, T;
            Historyresolution = Nstep, 
            zerofixpont       = false,
            affineinteration  = 1,
            Krylov_arg        = Krylov_arg)

        muaff, s0aff, solfix = affine(dpdploc; p = [ζ, ωn, k, w, K, Ωrel, RVA, RVF, f0, ϕshift])
        surface_pattern      = fnan.(g.(ϕ.(solfix.t, ϕshift)) .* real.(getindex.(solfix.u, [1])))
        surface_pattern      = filter(!isnan, surface_pattern)
        peak_to_peak         = maximum(surface_pattern) - minimum(surface_pattern)

        return abs(muaff[1][1]), peak_to_peak
    end

    Spek_aff[ii, jj], p2p[ii, jj] = fooTurningSSV(Ωrel, w)
    p2p_stab[ii, jj] = Spek_aff[ii, jj] < 1 ? p2p[ii, jj] : NaN
end
println("---------- End brute force ---------------")

# Turning speed at each angle
plot_phase()

# p2p(w,ϕ) surface at every points
plot(ϕshift_vec, w_vec*1e3, p2p'*1e6, st=:surface, c=reverse(cgrad(:hot)), zlim=(0,100), clims=(0, 10.0^(6 - 5)), xlabel="fi [rad]", ylabel="w [mm]", zlabel="peak to peak [um]")

# p2p(w,ϕ) surface at stable points
plot(ϕshift_vec, w_vec*1e3, p2p_stab'*1e6, st=:surface, c=reverse(cgrad(:hot)), zlim=(0,100), clims=(0, 10.0^(6 - 5)),title="3D Surface Plot", xlabel="fi [rad]", ylabel="w [mm]", zlabel="peak to peak [um]")

# p2p(w,ϕ) contour plot with heatmap at stable points
heatmap(ϕshift_vec, w_vec*1e3, p2p_stab'*1e6,  c=reverse(cgrad(:hot)), clims=(0, 10.0^(6 - 5)), colorbar=true, xlabel="fi [rad]", ylabel="w [mm]")
contour!(ϕshift_vec, w_vec*1e3, Spek_aff', levels=[1], lw=3, color=:gray)











# -------------------------------------------------------------------------------------------------------------
# ------------------------------------- p2p(Ωrel,w,ϕshift) isosurface -----------------------------------------
# -------------------------------------------------------------------------------------------------------------
RVF = 1
RVA = 0.1

ϕshift_vec  = range(0.0, 4*pi, 10)
w_vec       = range(0.01, 2.5, 10) / 1e3
Ωrel_vec    = range(1200.0, 2500.0, 10) * pi / 30

Aaff     = zeros(size(Ωrel_vec, 1), size(ϕshift_vec, 1), size(w_vec, 1))
Spek_aff = zeros(size(Ωrel_vec, 1), size(ϕshift_vec, 1), size(w_vec, 1))
p2p      = zeros(size(Ωrel_vec, 1), size(ϕshift_vec, 1), size(w_vec, 1))
p2p_stab = zeros(size(Ωrel_vec, 1), size(ϕshift_vec, 1), size(w_vec, 1))


println("---------- Starting brute force ---------------")
@time for (ii,Ωrel) in enumerate(Ωrel_vec), (jj,ϕshift) in enumerate(ϕshift_vec), (kk,w) in enumerate(w_vec)
    
    τmax  = 2pi / Ωrel / (1.0 - RVA) * 1.0
    Nstep = Int(ceil(τmax / dtapprox))

    function fooTurningSSV(Ωrel, w)

        dpdploc = dynamic_problemSampled(probTurningSSV, Solver_args, τmax, T;
            Historyresolution = Nstep, 
            zerofixpont       = false,
            affineinteration  = 1,
            Krylov_arg        = Krylov_arg)

        muaff, s0aff, solfix = affine(dpdploc; p=[ζ, ωn, k, w, K, Ωrel, RVA, RVF, f0, ϕshift])
        surface_pattern      = fnan.(g.(ϕ.(solfix.t, ϕshift)) .* real.(getindex.(solfix.u, [1])))
        surface_pattern      = filter(!isnan, surface_pattern)
        peak_to_peak         = maximum(surface_pattern) - minimum(surface_pattern)

        return abs(muaff[1][1]), peak_to_peak
    end

    Spek_aff[ii, jj, kk], p2p[ii, jj, kk] = fooTurningSSV(Ωrel, w)
    p2p_stab[ii, jj, kk] = Spek_aff[ii, jj, kk] < 1 ? p2p[ii, jj, kk] : NaN
end
println("---------- End brute force ---------------")

meshgrid    = Iterators.product(Ωrel_vec, ϕshift_vec, w_vec)
Ωrel_grid   = getindex.(meshgrid, 1)  
ϕshift_grid = getindex.(meshgrid, 2)
w_grid      = getindex.(meshgrid, 3)*1e3  
 
trace = isosurface(x=Ωrel_grid[:], y=ϕshift_grid[:], z=w_grid[:], value=Spek_aff[:], 
        opacity=0.9, isomin=1, isomax=1,surface_fill=1, surface_count=1, 
        colorscale=colors.YlOrRd_9, color=w_grid[:] )
layout = Layout(
    title="Isosurface",
    scene=attr(
        xaxis=attr(title="Ωrel [rad/s]"),    
        yaxis=attr(title="ϕshift [rad]"),  
        zaxis=attr(title="w [mm]")))
fig   = PlotlyJS.plot([trace], layout)