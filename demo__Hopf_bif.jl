#Demo: Delayed Nonline Oscill with nonlinearity
5 + 5

using Revise
using DDE_mapping

#using BenchmarkTools
using Plots
theme(:dark)#:vibrant:dracula:rose_pine
plotly()
#const PLOTS_DEFAULTS = Dict(:theme => :wong, :fontfamily => "Computer Modern", :label => nothing, :dpi => 600 )const PLOTS_DEFAULTS = Dict(:theme => :wong, :fontfamily => "Computer Modern", :label => nothing, :dpi => 600 )
default(size=(900, 900), titlefont=(15, "times"), legendfontsize=13, guidefont=(12, :white), tickfont=(12, :orange), guide="x", framestyle=:zerolines, yminorgrid=true, fontfamily="Computer Modern", label=nothing, dpi=600)

#using Profile
using StaticArrays
using DifferentialEquations

using LinearAlgebra

# using MDBM

using FileIO
using Suppressor

# using MDBM

using FileIO
# Governing equation

function DelayedNonlineOscill(u, h, p, t)
    # Parameters
    ζ, δ, b, τ, μ = p
    # Components of the delayed differential equation
    dx = u[2]
    ddx = -δ * u[1] - 2 * ζ * u[2] + b * h(p, t - τ)[1] - μ * u[2]^3 + μ * u[2]^5
    #ddx = -δ * u[1] - 2 * ζ * u[2] + b * h(p, t - τ)[1] - μ * u[2]^2*sign(u[2])
    # Update the derivative vector
    SA[dx, ddx]
end

Base.:+(a::SVector, b::Bool) = a .+ b
Base.:+(a::SVector, b::Float64) = a .+ b #TODO: where to put this?
ζ = 0.01         # damping coefficient
δ = 3.1#0.2          # nat. freq
b = -2.8#stable per.orbit
b = -3.15#unstable per.orbit
b = -1#Hopf starting...
b = 0.1#Hopf starting...
τ = 2pi#0.5#2pi          # Time delay
μ = 5.0
p = ζ, δ, b, τ, μ
#p = (ζ, ωn, k, τ,10.0)

# test simulation ---------------
#initial condition
u0 = SA[0.001, 0.0]
#history function
h(p, t) = SA[0.0; 0.0]

Tlongsim = 5000.2
Tend = 27.0
prob_long = DDEProblem{false}(DelayedNonlineOscill, u0, h, (0.0, Tlongsim), p; constant_lags=[τ])
#SciMLBase.AutoSpecialize
#SciMLBase.NoSpecialize
#SciMLBase.FullSpecialize

#Parameters for the solver as a Dict (it is necessary to collect it for later use)
Solver_args = Dict(:alg => MethodOfSteps(RK4()), :verbose => false, :reltol => 1e-7)#

@time sol = solve(remake(prob_long, p=(ζ, δ, b, τ, μ)); Solver_args...);
plot(sol)
##
#last period of the long simulation:
t_select_period = 0.0:0.01:Tend
t_select_delay = LinRange(0, τ, 200)
sol_period = sol(sol.t[end] .- t_select_period)
sol_delay = sol(sol.t[end] .- t_select_delay)

#plot the state phase of the last segment
plot(sol_period.t, sol_period[1, :])
plot!(sol_period.t, sol_period[2, :])
plot!(sol_delay.t, sol_delay[1, :])
plot!(sol_delay.t, sol_delay[2, :])
#plot the phase space (u - du)
plot(sol_delay[1, :], sol_delay[2, :])
plot!(sol_period[1, :], sol_period[2, :])


## ---------------- simulation max amplitude ----------------------
## parameters

bv = [LinRange(-2.0, 1.95, 100)..., LinRange(1.99, 2.05, 100)...]
norm_solperiod = similar(bv)
@time Threads.@threads for ib in eachindex(bv)
    println(ib)
    bloc = bv[ib]
    #prob_long = DDEProblem(DelayedNonlineOscill, u0, h, (0.0, Tlongsim), (ζ, δ, bloc, τ, μ); constant_lags=[τ])
    #sol = solve(prob_long; Solver_args...)#abstol,reltol

    sol = solve(remake(prob_long, p=(ζ, δ, bloc, τ, μ)); Solver_args...)#abstol,reltol

    #plot(sol)
    #last period of the long simulation:
    t_select_period = 0.0:0.01:Tend
    t_select_delay = eriod = 0.0:0.01:τ
    sol_period = sol(sol.t[end] .- t_select_period)
    sol_delay = sol(sol.t[end] .- t_select_delay)
    norm_solperiod[ib] = maximum(abs.(getindex.(sol_period.u, 1)))
end

plot(bv, norm_solperiod)
scatter!(bv, norm_solperiod)
##
#sol = solve(remake(prob_long, p=(ζ, δ, 2.010, τ, μ)); Solver_args...)#abstol,reltol
#plot(sol)





## ---------------- Affine mapping ---------------------------
#preparation for continuation:

using ForwardDiff
one_espilon_Dual = ForwardDiff.Dual{Float64}(0.0, 1.0)

using KrylovKit
Neig = 12#number of required eigen values
Krylov_arg = (Neig, :LM, KrylovKit.Arnoldi(tol=1e-25, krylovdim=3 + 20, verbosity=0, eager=true));

τmax = τ #maximal timedelay in the mapping
Nstep = 100 # discretization number of the mapping
Timeperiod = Tend # timeperiod of the mapping

#Creating the problem
dp_0_Tfix = dynamic_problemSampled(prob_long, Solver_args, τmax,
    Timeperiod; Historyresolution=Nstep,
    zerofixpont=true, affineinteration=0,
    Krylov_arg=Krylov_arg)

#bv_affine = LinRange(-2.0, 2.05, 201)#TODO: ez kell ha a balolbali szakaszt szeretném követeni
#bv_affine = LinRange(2.05, -2.0, 201)

bv_affine = [LinRange(-2.0, 1.95, 100)..., LinRange(1.99, 2.05, 100)...]
λ_μ₀ = Any[similar(bv_affine)...]

@time Threads.@threads for ib in eachindex(bv_affine)
    println(ib)
    bloc = bv_affine[ib]
    mu, saff, sol0 = affine(dp_0_Tfix; p=(ζ, δ, bloc, τ, μ))
    λ_μ₀[ib] = log.(abs.(mu[1])) / sol0.t[end]

end

plot(bv, norm_solperiod)
#scatter()
for k in 1:Neig
    plot!(bv_affine, getindex.(λ_μ₀, k))
end
plot!(legend=false)

## ---------------------------



##----------------------- DDE-Hopf amp stab -------------------


function condition(u, t, integrator) # Event when condition(u,t,integrator) == 0
    u[2]
end
function affect_short!(integrator)
    if integrator.t > 0.1
        terminate!(integrator)
    end
end
cb_short = ContinuousCallback(condition, affect_short!, nothing)

Solver_args_T_short = Dict(:alg => MethodOfSteps(RK4()), :verbose => false, :reltol => 1e-5,
    :callback => cb_short)#

using ForwardDiff
#TODO: ez nem hiszem, hogy jó megoldás
Base.:convert(::Type{Float64}, x::ForwardDiff.Dual{Float64,Float64,1}) = x.value

using KrylovKit
Neig = 8#number of required eigen values
Krylov_arg = (Neig, :LM, KrylovKit.Arnoldi(tol=1e-35, krylovdim=3 + 35, verbosity=0, eager=true));
#Creating the problem
Timeperiod = 20.0;#Tend



τmax = τ #maximal timedelay in the mapping
Nstep = 100 # discretization number of the mapping
Timeperiod = Tend # timeperiod of the mapping

#Creating the problem
dp_0_cb = dynamic_problemSampled(prob_long, Solver_args_T_short, τmax,
    Timeperiod; Historyresolution=Nstep,
    zerofixpont=false, affineinteration=3,
    Krylov_arg=Krylov_arg)





#------------------ initalization in a critical point --------------------------
bcrit = bv_affine[findfirst(x -> x < 0, getindex.(λ_μ₀, 1))]
b_H_start = 0.1;#-1.0#bcrit - 0.05;

mu_c, saff_c, sol0_c = affine(dp_0_Tfix; p=(ζ, δ, bcrit, τ, μ));
mu_c[1]
T = mu_c[3]
As = mu_c[2]
Acrit = mu_c[2][1]




plot(sol_delay[1, :], sol_delay[2, :])
plot!(sol_period[1, :], sol_period[2, :])
plot!(getindex.(Acrit, 1), getindex.(Acrit, 2), lw=1)
#plot!(dp_0_Tfix.StateSmaplingTime, getindex.(Acrit,1),lw=2)
#plot!(dp_0_Tfix.StateSmaplingTime, getindex.(Acrit,2),lw=2)

Acrtit_cb = LinMap(dp_0_cb, Acrit; p=(ζ, δ, b_H_start, τ, μ))[1]
plot!(getindex.(Acrtit_cb, 1), getindex.(Acrtit_cb, 2), lw=2)





mu_c, saff_c, sol0_c = affine(dp_0_cb, Acrtit_cb * 10.0; p=(ζ, δ, b_H_start, τ, μ));
scatter!(getindex.(saff_c, 1), getindex.(saff_c, 2), lw=2)
scatter!(sol0_c[1, :], sol0_c[2, :])


# ----------- brute force naive continuation ----------------

#ustart = rand(typeof(dp_0_cb.Problem.u0), Nstep)
#ustart = ustart .* 0.0 .+ 1.5
ustart = saff_c
#bv_affine_H = LinRange(-2.0, -0.951233, 27)
#bv_affine_H = LinRange(-1.0, -4.0, 37)
#bv_affine_H = LinRange(1.2, 0.08, 27)
bv_affine_H = LinRange(b_H_start, -2.0, 150)
bv_affine_H = [LinRange(b_H_start, -1.8, 10)..., LinRange(-1.81, -1.96, 40)..., LinRange(-1.96, -1.98, 100)...]

bv_affine_H = [LinRange(b_H_start, 2, 100)...]
bv_affine_H = [b_H_start:0.01:2.0...]

bv_affine_H = [LinRange(b_H_start, 1.95, 30)..., LinRange(1.96, 2.01, 30)...]

λ_μ₀_Hopf = Any[similar(bv_affine_H)...]
Amp_H = Any[similar(bv_affine_H)...]
T_period_H = Any[similar(bv_affine_H)...]
#
#plot()
#Threads.@threads

@suppress_err begin
    @time for ib in eachindex(bv_affine_H)
        println(ib)
        @show bloc = deepcopy(bv_affine_H[ib])
        mu, saff, sol0 = affine(dp_0_cb, ustart; p=(ζ, δ, bloc, τ, μ))
        #mu, saff, sol0 = affine(dp_Hopf_callback; p=(ζ, δ, bloc, τ, μ))
        λ_μ₀_Hopf[ib] = log.(abs.(mu[1])) / sol0.t[end]
        #Amp[ib] = saff
        T_period_H[ib] = sol0.t[end]
        Amp_H[ib] = maximum(abs.(getindex.(saff, 1)))

        ustart = saff #TODO: enélkül nem áll be szépen (de csak azért nem csinálta, mert affineinteration nem érvényesült)

        #   plot!(sol0[1, :], sol0[2, :])
        #   aaa=plot!(getindex.(saff, 1), getindex.(saff, 2), marker=:circle, markersize=1, lab="")#
        #   display(aaa)


    end
end
#plot!(legend=false)

##--------
plot(bv, norm_solperiod)
marcolor = sign.(getindex.(λ_μ₀_Hopf, 1))
#, color=:bamako
scatter!(bv_affine_H, Amp_H, zcolor=marcolor, lw=0, markerstrokewidth=0)
plot!(legend=false)
plot!(ylim=(-0.3, 1000.0))
plot!(ylim=())



for k in 1:Neig
    plot!(bv_affine_H, getindex.(λ_μ₀_Hopf, k))
    #for ib in eachindex(bv_affine_H)
    # scatter!([bv_affine_H[ib]], [λ_μ₀_Hopf[ib][k]])
    #end
end
plot!()
plot!(ylim=(-0.5, 1.0))

fig_bif = plot!(bv_affine_H, T_period_H ./ 10.0, lw=3)


# ------------------------------------------------------------------------------
# TODO: csak itt tartok -  tesztelgetés
fig_normerror = scatter(yaxis=:log)
fig_abs_mus = scatter(yaxis=:log)
fig_path = plot()
# Pseudo-archlength method
5 + 5
b_H_start = 0.1
λ_start = b_H_start
p = (ζ, δ, b_H_start, τ, μ)
p_start = p
#initialization of a step:


mu_c, s0_start, sol0_c = affine(dp_0_cb, Acrtit_cb * 10.0; p=p_start);
mu_c, s0_start, sol0_c = affine(dp_0_cb, s0_start; p=p_start);


Δλ = 0.05#05
#Δλ = 0.005
#Δλ = 0.05
p = p .+ (0.0, 0.0, Δλ, 0.0, 0.0)
mu_c, s0, sol0_c = affine(dp_0_cb, s0_start; p=p);
mu_c, s0, sol0_c = affine(dp_0_cb, s0; p=p);



scatter!(fig_bif, [p_start[3]], [maximum(abs.(getindex.(s0_start, 1)))], markersize=2)
scatter!(fig_bif, [p[3]], [maximum(abs.(getindex.(s0, 1)))], markersize=2)
#scatter!([p[3]], [maximum(abs.(getindex.(s0, 1)))])
# # #mu_c_m1, saff_c_m1, sol0_c_m1 = affine(dp_0_cb, s0_start; p=(ζ, δ, b_H_start - 0.01, τ, μ));
# # #scatter(getindex.(saff_c, 1), getindex.(saff_c, 2), lw=2)
# # #scatter!(sol0_c[1, :], sol0_c[2, :])
# # #
# # #maximum(abs.(getindex.(saff_c, 1)))
# # #maximum(abs.(getindex.(saff_c_m1, 1)))
# # #Δu_Δλ = (saff_c .- saff_c_m1) ./ 0.01
# # #
# # #
# # #f(x) = affine(dp_0_cb, Acrtit_cb * 10.0; p=(ζ, δ, b_H_start + x, τ, μ))[2]
# # #df(x) = (f(x) - f(x - 0.01)) ./ 0.01
# # #ddd = df(0.00)
# # #
# # #
# # #
# # #Δλ = 0.03
# # #p = p .+ (0.0, 0.0, Δλ, 0.0, 0.0)
# # #s0 = s0_start .+ Δu_Δλ
#λscale=10.0
#for _ in 1:5#32

#using JLD2
#@save "bifurcation_till_Fold_point.jld2" 
#@load "bifurcation_till_Fold_point.jld2"
norm_limit = 1e-5
Niteration = 4
iteration_gola=7
maxΔλ=0.2


@suppress_err begin
    λscale = 10.0
    lamscale = 1.0
    Nconti = 1
    for Nconti in 1:(5*340)#34
        println("---------------------------")
        @show Nconti
        Δu = s0 - s0_start
        @show Δλ = p[3] - λ_start

        s0_start = deepcopy(s0)
        p_start = deepcopy(p)
        @show λ_start = p[3]
        #scatter!([p_start[3]], [maximum(abs.(getindex.(s0_start, 1)))],markersize=5)
        @show Niteration
        @show lamscale = 1.3 .^ (iteration_gola .- Niteration)
        lamscale=minimum([lamscale,1.5])
        # # #  if p[3] > 2.0
        # # #      break
        # # #  end
        # # #  if p[3] > 1.85
        # # #      @show lamscale = 1.0#0.7#0.0#0.7
        # # #  else
        # # #      @show lamscale = 1.0
        # # #  end
        @show lamscale=(Δλ * lamscale)>maxΔλ ? 1.0 : lamscale

        s0 .+= Δu * lamscale
        #s0 = s0 + Δu_Δλ * Δλ
        p = p .+ (0.0, 0.0, Δλ * lamscale, 0.0, 0.0)
        #mu, saff, sol0 = affine(dp_0_Tfix; p=(ζ, δ, bloc+one_espilon_Dual, τ, μ))
        dp = dp_0_cb
        Nstep = size(dp.StateSmaplingTime, 1)
        #s0 = rand(typeof(dp.Problem.u0), Nstep)
        #s0 .*= 0.0
        #s0 = Acrtit_cb * 10.0
        @warn "perturbation to see the convergence"
        #s0 .+= 0.001 
        #s0 .*= 1.05

        scatter!(fig_bif, [p[3]], [maximum(abs.(getindex.(s0, 1)))], markersize=3)
        scatter!(fig_bif, xlim=(0.09, 0.25), ylim=(0.05, 0.12))



        # iteration within a step -------
        Niteration = 0
        eigval=Any[]
        for kk_out in 1:5
            #pDual=p .* 0.0 .+ ForwardDiff.Dual{Float64}(0.0, 1.0)
            pDual = p .+ (0.0, 0.0, ForwardDiff.Dual{Float64}(0.0, 1.0), 0.0, 0.0)

            v0 = LinMap(dp, s0; p=p)[1]
            v0_dual = LinMap(dp, s0; p=pDual)[1]
            v0 = valuepart.(v0_dual)
            dv0dλ = partialpart.(v0_dual)


            #v02 = LinMap(dp, s0; p=(ζ, δ, b_H_start, τ, μ))[1]
            #v01 = LinMap(dp, s0; p=(ζ, δ, b_H_start - 0.001, τ, μ))[1]
            #dvdl = (v02 - v01) / 0.001
            #dv0dλ - dvdl

            norm(s0 - v0)

            scatter!(fig_normerror, [p[3]], [norm(s0 - v0)], markersize=3, yaxis=:log)
            TheMapping(s) = partialpart.(LinMap(dp, s * one_espilon_Dual + s0; p=p)[1] - v0)
            Nstep = size(dp.StateSmaplingTime, 1)
            s_start = rand(typeof(dp.Problem.u0), Nstep)

            mus = getindex(schursolve(TheMapping, s_start, dp.Krylov_arg...), [3, 2, 1])


            eigval = mus[1]
            eigvec = mus[2]


            v0 = LinMap(dp, s0; p=p)[1]

            #scatter!([p[3]], [maximum(abs.(getindex.(s0, 1)))])
            scatter!(fig_bif, xlim=(0.09, 0.25), ylim=(0.05, 0.12))

            #a0 = real.(find_fix_pont(s0, v0, mus[1], mus[2]))
            ## println("find_fix_pont_end")
            ## s0 = find_fix_pont(s0, LinMap(dp, s0; p=p), mus[1], mus[2])
            #@show normerror = norm(s0 - v0)
            #s0 = a0;
            #scatter!([p[3]],[maximum(abs.(getindex.(s0, 1)))],xlim=(0.09,0.13),ylim=(0.05,0.07))

            begin
                for kkkk in 1:5
                    
                Niteration += 1
                    #a0 = real.(find_fix_pont(s0, v0, mus[1], mus[2],dv0dλ,Δu_Δλ))
                    x = (v0 - s0)

                    scatter!(fig_normerror, [p[3]], [norm(s0 - v0)], markersize=1, yaxis=:log)
                    # println("------------------------------------")
                    # println(AtA)
                    Atx = [dot(eigvec[i], x) for i in 1:size(eigvec, 1)]



                    #ci = Atx #TODO: ez ugyan azt adja Schur esetén!!!
                    #ci_mu = (ci .* ((eigval) ./ (eigval .- 1.0)))#TODO: Szabad ezt csinálni, a Schur-nál, nem a sajátértékkel kellenen skálázni... (vagy az pont kiesik valós függvényeknél???)
                    #fix_v = v0 - mapreduce(x -> x[1] * x[2], +, zip(eigvec, ci_mu))
                    #Δλ = 0.0
                    #error("itt tartok, valahol itt kell előjeleket kereseni talán")

                    Atdvdlam = [dot(eigvec[i], dv0dλ) for i in 1:size(eigvec, 1)]
                    #Δu_ΔλAt = [dot(Δu_Δλ, eigvec[i]) for i in 1:size(eigvec, 1)]
                    #T_jac = vcat(hcat(diagm(eigval .- 1.0), -Atdvdlam),
                    #    hcat(-Δu_ΔλAt', 100.0))
                    ΔuAt = [dot(Δu, eigvec[i]) for i in 1:size(eigvec, 1)]
                    T_jac = vcat(hcat(diagm(eigval .- 1.0), -Atdvdlam),
                        hcat(-ΔuAt', λscale * Δλ))
                    ci_arch = T_jac \ vcat(Atx, 0.0)
                    ci = ci_arch[1:end-1]
                    Δλ_loc = real.(ci_arch[end])
                    ci_mu_ach = ci .* (eigval)
                    ci_mu = ci_mu_ach
                    #A=transpose(mapreduce(permutedims, vcat, eigvec))
                    #fix_v = v0 - A * (((A'A) \ (A' * x)) .* ((eigval) ./ (eigval .- 1.0)))

                    fix_v = v0 - mapreduce(x -> x[1] * x[2], +, zip(eigvec, ci_mu))
                    #fix_v = s0 - mapreduce(x -> x[1] * x[2], +, zip(eigvec, ci))
                    s0 = real.(fix_v)

                    v0_dual = LinMap(dp, s0; p=pDual)[1]
                    v0 = valuepart.(v0_dual)
                    dv0dλ = partialpart.(v0_dual)


                    p = p .+ (0.0, 0.0, Δλ_loc, 0.0, 0.0)

                    v0 = LinMap(dp, s0; p=p)[1]
                    if norm(s0 - v0) < norm_limit
                        break
                    end

                    scatter!(fig_bif, [p[3]], [maximum(abs.(getindex.(s0, 1)))], markersize=1)
                    scatter!(fig_bif, xlim=(2.002, 2.02), ylim=(0.43, 0.45))
                end
                if norm(s0 - v0) < norm_limit
                    break
                end
            end
            scatter!(fig_bif, [p[3]], [maximum(abs.(getindex.(s0, 1)))], markersize=2)
            #scatter!(xlim=(0.09,0.13),ylim=(0.05,0.07))

        end
        scatter!(fig_bif, [p[3]], [maximum(abs.(getindex.(s0, 1)))], m=:cross, markersize=4)
        #scatter!(xlim=(0.09,0.13),ylim=(0.05,0.07))
        #aaa = scatter!(fig_bif,xlim=(0.00, 2.3), ylim=(0.00, 0.8))
        scatter!(fig_abs_mus, p[3] .* ones(size(eigval)), abs.(eigval), markersize=3, m=:circle, yaxis=:log,xlabel="λ",ylabel="asb_mu")
        
        
        v0, sol = LinMap(dp, s0; p=p)
        plot!(fig_path,getindex.(v0, 1), getindex.(v0, 2))#,xlim=(),ylim=()
        plot!(fig_path,sol[1, :], sol[2, :],xlabel="u1",ylabel="u2",title="path", lw=2)

        scatter!(fig_normerror, [p[3]], [norm(s0 - v0)], markersize=5, m=:cross, yaxis=:log,xlabel="λ",ylabel="Norm error")
        scatter!(fig_bif, xlim=(0.00, 2.3), ylim=(0.00, 0.8),xlabel="λ",ylabel="norm(u)")
        out = plot(fig_bif, fig_normerror, fig_abs_mus,fig_path,layout=(2, 2))
        display(out)
    end
end
scatter!(fig_bif, xlim=(0.00, 2.3), ylim=(0.00, 0.8))
scatter!(fig_normerror, xlim=(0.00, 2.3))
plot(fig_bif, fig_normerror, layout=(2, 1))

# using JLD2
# #@save "bifurcation_till_Fold_point.jld2" 
# @load "bifurcation_till_Fold_point.jld2"
































#### 
#### # ----------------------- creating stability chart -------------------------
####
#### Krylov_arg = (1, :LM, KrylovKit.Arnoldi(tol=1e-12, krylovdim=8 + 5, verbosity=0));
####
#### dpMathieu = dynamic_problemSampled(probMathieu, Solver_args, τmax,
####     Timeperiod; Historyresolution=Nstep,
####     zerofixpont=false, affineinteration=1,
####     Krylov_arg=Krylov_arg)
####
#### ## ------------------ Peak-2-Peak ampNorm color, and spectral radius - BRUTE FORCE in the plane of δ-ϵ  ----------------
####
#### b = 0.05 # delay parameter, Note: b=0 provide the traditional stability chart of the Mathieu equations
#### δv = -1.0:0.051:5.0 # initial grid in x direction
#### ϵv = -0.001:0.05:10.0 # initial grid in y direction
#### Aaff = zeros(size(ϵv, 1), size(δv, 1))
#### Spek_aff = zeros(size(ϵv, 1), size(δv, 1))
####
#### @time Threads.@threads for j in 1:size(δv, 1)
####     @inbounds δ = δv[j]
####     Threads.@threads for i in 1:size(ϵv, 1)
####         @inbounds ϵ = ϵv[i]
####         muaff, s0aff, sol0 = affine(dpMathieu; p=(ζ, δ, ϵ, b, τ, T))
####         Aaff[i, j] = norm(getindex.(s0aff, 1)) # norm of the motion
####         # Aaff[i,j]= maximum(abs.(getindex.(s0aff,1))) #maximum amplitude of the orbit
####         Spek_aff[i, j] = maximum(abs.(muaff[1]))
####     end
#### end
####
#### #Plotting the maximal amplitud on the stable domain only
#### Aaffsat = deepcopy(Aaff);
#### Aaffsat[Spek_aff.>1.0] .= 0.0;#eliminate the positions of instable case
#### heatmap(δv, ϵv, log.(Aaffsat))
####
####
####
#### #Plotting the maximal amplitud on the stable domain only
#### Spek_affsat = deepcopy(Spek_aff);
#### Spek_affsat[Spek_affsat.>1.0] .= 0.0;
#### heatmap(δv, ϵv, log.(Spek_affsat))
####
####
#### #------------------ Stability boundary - MDBM -----------------
####
#### b = 0.05;
#### ax1 = Axis(0:1:5, "δ") # initial grid in x direction
#### ax2 = Axis(-0.001:2.0:10.0, "ϵ") # initial grid in y direction
#### function fooDelay(δ, ϵ)
####     ABSmuMax = spectralradius(dpMathieu; p=(ζ, δ, ϵ, b, τ, T))
####     #return ABSmuMax - 1.0
####     return log(ABSmuMax)
#### end
#### mymdbm = MDBM_Problem(fooDelay, [ax1, ax2])
#### @time MDBM.solve!(mymdbm, 5)
#### #points where the function foo was evaluated
#### x_eval, y_eval = getevaluatedpoints(mymdbm)
#### x_sol, y_sol = getinterpolatedsolution(mymdbm)
#### #scatter(x_eval,y_eval,markersize=1)
#### #scatter!(x_sol,y_sol,markersize=2)
#### scatter!(x_sol, y_sol, markersize=1)
####
####
####
#### #--------------------------
####
####
####
#### ## ------------------ Peak-2-Peak ampNorm color, and spectral radius - BRUTE FORCE  in the plane of δ-b-----------------
####
#### Krylov_arg = (1, :LM, KrylovKit.Arnoldi(tol=1e-12, krylovdim=8 + 5, verbosity=0));
####
#### dpMathieu = dynamic_problemSampled(probMathieu, Solver_args, τmax,
####     Timeperiod; Historyresolution=Nstep,
####     zerofixpont=false, affineinteration=1,
####     Krylov_arg=Krylov_arg)
####
####
#### δv = 0:0.051:10 # initial grid in x direction
#### bv = -1.501:0.05:1.5 # initial grid in y direction
#### Aaff = zeros(size(bv, 1), size(δv, 1))
#### Spek_aff = zeros(size(bv, 1), size(δv, 1))
####
#### @time Threads.@threads for j in 1:size(δv, 1)
####     @inbounds δ = δv[j]
####     Threads.@threads for i in 1:size(bv, 1)
####         @inbounds b = bv[i]
####         muaff, s0aff, sol0 = affine(dpMathieu; p=(ζ, δ, ϵ, b, τ, T))
####         Aaff[i, j] = norm(getindex.(s0aff, 1)) # norm of the motion
####         # Aaff[i,j]= maximum(abs.(getindex.(s0aff,1))) #maximum amplitude of the orbit
####         Spek_aff[i, j] = maximum(abs.(muaff[1]))
####     end
#### end
####
#### #Plotting the maximal amplitud on the stable domain only
#### Aaffsat = deepcopy(Aaff);
#### Aaffsat[Spek_aff.>1.0] .= 0.0;#eliminate the positions of instable case
#### heatmap(δv, bv, log.(Aaffsat))
####
####
####
#### #Plotting the maximal amplitud on the stable domain only
#### Spek_affsat = deepcopy(Spek_aff);
#### Spek_affsat[Spek_affsat.>1.0] .= 0.0;
#### heatmap(δv, bv, log.(Spek_affsat))
####
####
#### #------------------ Stability boundary - MDBM -----------------
####
####
#### ax1 = Axis(0:2:10, "δ") # initial grid in x direction
#### ax2 = Axis(-1.5:1.4:1.5, "b") # initial grid in y direction
#### function fooDelay(δ, b)
####     ABSmuMax = spectralradius(dpMathieu; p=(ζ, δ, ϵ, b, τ, T))
####     return ABSmuMax - 1.0
#### end
#### mymdbm = MDBM_Problem(fooDelay, [ax1, ax2])
#### @time MDBM.solve!(mymdbm, 4)
#### #points where the function foo was evaluated
#### x_eval, y_eval = getevaluatedpoints(mymdbm)
#### x_sol, y_sol = getinterpolatedsolution(mymdbm)
#### #scatter(x_eval,y_eval,markersize=1)
#### #scatter!(x_sol,y_sol,markersize=2)
#### scatter!(x_sol, y_sol, markersize=1)
####