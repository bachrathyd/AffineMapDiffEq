#TODO_ valamiért nem megy, szétkenődik a megoldás, de lehet, hogy ez így van rendben!??!?!
#Bending vibrations

5 + 5


using ModelingToolkit, MethodOfLines, OrdinaryDiffEq, DifferentialEquations, DomainSets

using Plots

# Define constants (replace with your specific values)
E_modulus = 210e9# Steal
a_1 = 20e-3#size of crosssection horizontal
b_1 = 3e-3#size of crosssection vertivcal
ρ_densitiy = 7800#
Area = a_1 * b_1
Inertia = a_1 * b_1^3 / 12

EI = E_modulus * Inertia
ρA = Area * ρ_densitiy



η0 = 0#.000000001
η(x) = η0#damping

L = 10.0 # Length of the beam
Tend = 1.5 # end time

@parameters x t
@variables w(..)
Dt = Differential(t)
Dtt = Differential(t)^2
Dx = Differential(x)
Dxx = Differential(x)^2
Dxxx = Differential(x)^3


ω = 10.0*2*pi
q(x, t) = 0.0 * cos(ω * t)

#w0(x, t) = 0.01 * x^2
w0(x) = 0.01 * exp(-((x - L / 3) / L * 150/4)^2)
w0(x, t) = 1.0*w0(x)
BC_time(t) = 0.01 * exp(-((t -0.02) / 0.005)^2)
#BC_time(x, t)  = BT_time(t)
#BC_time(t) = 0.01 * sin(ω*t)
BC_time(x, t)  = BC_time(t)
#display(plot(w0.(LinRange(0, L, 500))))
#display(plot(BC_time.(LinRange(0, Tend, 500))))


#eq = [ρA * Dtt(w(x, t)) + η(x)*Dt(w(x, t)) + Dxx(EI*(Dxx(w(x, t))))~ + q(x, t) ]
eq = [ρA * Dtt(w(x, t)) + η(x) * Dt(w(x, t)) + EI * Dxx((Dxx(w(x, t)))) ~ +q(x, t)]




# Define spatial domain and grid

x_min = t_min = 0.0
x_max = L
t_max = Tend

domains = [x ∈ Interval(x_min, x_max),
    t ∈ Interval(t_min, t_max)]


#if false#true
# #Boundary conditions #Clamped-Free
#bcs = [ w(x, 0) ~ w0(x, t),
#        Dt(w(x, 0)) ~ 0.0,
#        w(0, t) ~  0.0,        
#        Dx(w(0, t)) ~  0.0,    
#        Dxx(w(x_max, t)) ~ 0.0,    
#        Dxxx(w(x_max, t)) ~ 0.0] 
#else
# # Boundary conditions - #Clamped-Clamped
bcs = [w(x, 0) ~ w0(x, t),
    Dt(w(x, 0)) ~ 0.0,
    w(0, t) ~ 0.0,#BC_time(t),
    Dx(w(0, t)) ~ 0.0,
    w(x_max, t) ~ 0.0,
    Dx(w(x_max, t)) ~ 0.0]
#end




@named pdesys = PDESystem(eq, bcs, domains, [x, t], [w(x, t)])


N = 120  # Number of spatial grid points
order = 4  # Order of the spatial discretization
discretization = MOLFiniteDifference([x => N], t, approx_order=order, jac = true, sparse = true)


println("Discretization:")
@time prob = discretize(pdesys, discretization)



println("Solve:")
@time sol = solve(prob, Tsit5(), reltol=1e-3, dtmax=1e-2)
#@time sol = solve(prob, TRBDF2(), reltol=1e-3, dtmax=1e-2)

discrete_x = sol[x]
discrete_t = sol[t]

solw = sol[w(x, t)]
@show size(solw)

#plotly()
gr()
display(heatmap(discrete_t[1:100:end],discrete_x,solw[:,1:100:end]))
#display(heatmap(discrete_t[1:10:end],discrete_x,solw[:,1:10:end]))
display(heatmap(discrete_t[1:1:end],discrete_x,solw[:,1:1:end]))

##

# Make an animation
anim = @animate for k in 1:10:(length(discrete_t)÷1)
    plot(discrete_x,solw[1:end, k], title="$(discrete_t[k])",ylim=(-0.01,0.01)) # 2:end since end = 1, periodic condition
end
display(gif(anim, "Beam_vibration_1.gif", fps = 25))

println("-------------------------------- DONE --------------------------------")

##
@time sol_no_wrap = solve(prob, Tsit5(), reltol=1e-6, dtmax=1e-3; wrap = Val(false))
plot(sol_no_wrap(0.1)[501:end])
anim = @animate for tloc in LinRange(0, Tend, 500)
    plot(sol_no_wrap(tloc)[500:end], title="$tloc",ylim=(-0.01,0.01)) # 2:end since end = 1, periodic condition
end
display(gif(anim, "Beam_vibration_time.gif", fps = 25))
##


#TODO:  Eddig kellene működnie, de kell majd velemi ellnőrzés.












































5 + 5
using Interpolations
using StaticArrays
using DifferentialEquations
using Plots
plotly()
#gr()
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


using Revise
using DDE_mapping




Base.:+(a::SVector, b::Bool) = a .+ b

##
# -------------------------------------------------------------------------------------------------------------
# ------------------------------------- Testing affine mapping ------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
T=0.01
    Krylov_arg = (5, :LM, KrylovKit.Arnoldi(tol=1e-15, krylovdim=10, verbosity=3))
    Solver_args = Dict(:alg =>MethodOfSteps(Tsit5()), :verbose => false, 
    :reltol => 1e-8, :maxiters => 1e7, :saveat => LinRange(0, T, 250),wrap => Val(false))

    τmax = 0.001
    dtapprox = τmax / 10

    @time sol_no_wrap = solve(prob, Tsit5(), reltol=1e-6, dtmax=1e-3; wrap = Val(false))
    plot(sol_no_wrap(0.1)[501:end])
    anim = @animate for tloc in LinRange(0, Tend, 500)
        plot(sol_no_wrap(tloc)[500:end], title="$tloc",ylim=(-0.01,0.01)) # 2:end since end = 1, periodic condition
    end
    display(gif(anim, "Beam_vibration_time.gif", fps = 25))



    @show Nstep = Int(ceil(τmax / dtapprox)) 

    probTurningSSV = DDEProblem(prob, u0, h, [0, T], p)
    dpdp = dynamic_problemSampled(probTurningSSV, Solver_args, τmax, T;
        Historyresolution=Nstep,
        zerofixpont=false,
        affineinteration=2,
        Krylov_arg=Krylov_arg)

    @time muaff, s0aff, solfix = affine(dpdp; p=[ζ, ωn, k, w, K, Ωrel, RVA, RVF, f0, ϕshift])
    println(abs.(muaff[1]))
    fnan(y) = y == 0 ? NaN : y
    plot!(solfix.t, real.(getindex.(solfix.u, 1)))
    plot!(solfix.t, fnan.(g.(ϕ.(solfix.t, ϕshift)) .* real.(getindex.(solfix.u, 1))), lw=3, ls=:dot)

plot!()

##


# -------------------------------------------------------------------------------------------------------------
# ------------------------------------- p2p(w,ϕshift) surface and contour plot with heatmap -------------------
# -------------------------------------------------------------------------------------------------------------
plot()
for RVA in 10.0 .^ (-5:0.1:-1)
    # begin
    #    RVA = 0.1
    println(RVA)
    RVF = 1
    Krylov_arg = (2, :LM, KrylovKit.Arnoldi(tol=1e-10, krylovdim=5, verbosity=0))
    Solver_args = Dict(:alg => MethodOfSteps(RK4()), :verbose => false, :reltol => 1e-5, :maxiters => 1e7,)

    Ωrel_vec = 230

    ϕshift_vec = [1.5]#,1.8]
    #ϕshift_vec  = range(0.0, 4*pi,10)
    w_vec = range(0.001, 2.5, 50) / 1e3
    #w_vec = range(0.0, 1.5, 15) / 1e3

    Aaff = zeros(size(ϕshift_vec, 1), size(w_vec, 1))
    Spek_aff = zeros(size(ϕshift_vec, 1), size(w_vec, 1))
    p2p = zeros(size(ϕshift_vec, 1), size(w_vec, 1))
    p2p_stab = zeros(size(ϕshift_vec, 1), size(w_vec, 1))


    println("---------- Starting brute force ---------------")
    @time for (ii, ϕshift) in enumerate(ϕshift_vec), (jj, w) in enumerate(w_vec)

        T = maximum([2pi / (Ωrel * RVF), 2pi / Ωrel])
        τmax = 2pi / Ωrel / (1.0 - RVA) * 1.0
        Nstep = Int(ceil(τmax / dtapprox)) * 5

        function fooTurningSSV(Ωrel, w)

            dpdploc = dynamic_problemSampled(probTurningSSV, Solver_args, τmax, T;
                Historyresolution=Nstep,
                zerofixpont=false,
                affineinteration=2,
                Krylov_arg=Krylov_arg)

            muaff, s0aff, solfix = affine(dpdploc; p=[ζ, ωn, k, w, K, Ωrel, RVA, RVF, f0, ϕshift])
            surface_pattern = fnan.(g.(ϕ.(solfix.t, ϕshift)) .* real.(getindex.(solfix.u, [1])))
            surface_pattern = filter(!isnan, surface_pattern)
            peak_to_peak = maximum(surface_pattern) - minimum(surface_pattern)

            return abs(muaff[1][1]), peak_to_peak
        end

        Spek_aff[ii, jj], p2p[ii, jj] = fooTurningSSV(Ωrel, w)
        p2p_stab[ii, jj] = Spek_aff[ii, jj] < 1 ? p2p[ii, jj] : NaN
    end
    println("---------- End brute force ---------------")

    # Turning speed at each angle
    #plot_phase()
    #p2p[p2p .> 0.08] .= 0.0
    # p2p(w,ϕ) surface at every points
    #plot(ϕshift_vec, w_vec*1e3, p2p'*1e6, st=:surface, c=reverse(cgrad(:hot)), zlim=(0,100), clims=(0, 10.0^(6 - 5)), xlabel="fi [rad]", ylabel="w [mm]", zlabel="peak to peak [um]")

    #plot(w_vec*1e3, (p2p[:]*1e6), zlim=(0,100), xlabel="w [mm]", ylabel="peak to peak [um]")
    plot!(w_vec * 1e3, (p2p[:] * 1e5), xlabel="w [mm]", ylabel="peak to peak [um]", ylim=(0, 20))

    aa = plot!(w_vec * 1e3, Spek_aff[:], xlabel="w [mm]", ylabel="stab_mu", ylim=(0, 20), lw=2, ls=:dot)
    display(aa)
    # p2p(w,ϕ) surface at stable points
    #plot(ϕshift_vec, w_vec*1e3, p2p_stab'*1e6, st=:surface, c=reverse(cgrad(:hot)), zlim=(0,100), clims=(0, 10.0^(6 - 5)),title="3D Surface Plot", xlabel="fi [rad]", ylabel="w [mm]", zlabel="peak to peak [um]")

    ## p2p(w,ϕ) contour plot with heatmap at stable points
    #heatmap(ϕshift_vec, w_vec*1e3, p2p_stab'*1e6,  c=reverse(cgrad(:hot)), clims=(0, 10.0^(6 - 5)), colorbar=true, xlabel="fi [rad]", ylabel="w [mm]")
    #contour!(ϕshift_vec, w_vec*1e3, Spek_aff', levels=[1], lw=3, color=:gray)
end
plot!(ylim=(0, 20))



##

# -------------------------------------------------------------------------------------------------------------
# ------------------------------------- p2p(Ωrel,w,ϕshift) isosurface -----------------------------------------
# -------------------------------------------------------------------------------------------------------------
RVF = 1
RVA = 0.1

ϕshift_vec = range(0.0, 4 * pi, 10)
w_vec = range(0.01, 2.5, 10) / 1e3
Ωrel_vec = range(1200.0, 2500.0, 10) * pi / 30

Aaff = zeros(size(Ωrel_vec, 1), size(ϕshift_vec, 1), size(w_vec, 1))
Spek_aff = zeros(size(Ωrel_vec, 1), size(ϕshift_vec, 1), size(w_vec, 1))
p2p = zeros(size(Ωrel_vec, 1), size(ϕshift_vec, 1), size(w_vec, 1))
p2p_stab = zeros(size(Ωrel_vec, 1), size(ϕshift_vec, 1), size(w_vec, 1))


println("---------- Starting brute force ---------------")
@time for (ii, Ωrel) in enumerate(Ωrel_vec), (jj, ϕshift) in enumerate(ϕshift_vec), (kk, w) in enumerate(w_vec)

    τmax = 2pi / Ωrel / (1.0 - RVA) * 1.0
    Nstep = Int(ceil(τmax / dtapprox))

    function fooTurningSSV(Ωrel, w)

        dpdploc = dynamic_problemSampled(probTurningSSV, Solver_args, τmax, T;
            Historyresolution=Nstep,
            zerofixpont=false,
            affineinteration=1,
            Krylov_arg=Krylov_arg)

        muaff, s0aff, solfix = affine(dpdploc; p=[ζ, ωn, k, w, K, Ωrel, RVA, RVF, f0, ϕshift])
        surface_pattern = fnan.(g.(ϕ.(solfix.t, ϕshift)) .* real.(getindex.(solfix.u, [1])))
        surface_pattern = filter(!isnan, surface_pattern)
        peak_to_peak = maximum(surface_pattern) - minimum(surface_pattern)

        return abs(muaff[1][1]), peak_to_peak
    end

    Spek_aff[ii, jj, kk], p2p[ii, jj, kk] = fooTurningSSV(Ωrel, w)
    p2p_stab[ii, jj, kk] = Spek_aff[ii, jj, kk] < 1 ? p2p[ii, jj, kk] : NaN
end
println("---------- End brute force ---------------")

meshgrid = Iterators.product(Ωrel_vec, ϕshift_vec, w_vec)
Ωrel_grid = getindex.(meshgrid, 1)
ϕshift_grid = getindex.(meshgrid, 2)
w_grid = getindex.(meshgrid, 3) * 1e3

trace = isosurface(x=Ωrel_grid[:], y=ϕshift_grid[:], z=w_grid[:], value=Spek_aff[:],
    opacity=0.9, isomin=1, isomax=1, surface_fill=1, surface_count=1,
    colorscale=colors.YlOrRd_9, color=w_grid[:])
layout = Layout(
    title="Isosurface",
    scene=attr(
        xaxis=attr(title="Ωrel [rad/s]"),
        yaxis=attr(title="ϕshift [rad]"),
        zaxis=attr(title="w [mm]")))
fig = PlotlyJS.plot([trace], layout)