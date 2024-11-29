#Demo: Delayed Nonline Oscill with nonlinearity
5 + 5

using Revise
using DDE_mapping

using BenchmarkTools
using Plots
theme(:dark)#:vibrant:dracula:rose_pine

plotly()
using Profile
using StaticArrays
using DifferentialEquations

using LinearAlgebra

using MDBM


# Governing equation

function DelayedNonlineOscill(u, h, p, t)
    # Parameters
    ζ, δ, b, τ, μ = p
    # Components of the delayed differential equation
    dx = u[2]
    ddx = -δ * u[1] - 2 * ζ * u[2] + b * h(p, t - τ)[1] - μ * u[2]^3+ μ * u[2]^5
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

bv = LinRange(-2.0, 2.05, 110)
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
##





## ---------------- Affine mapping ---------------------------
#preparation for continuation:

using ForwardDiff
one_espilon_Dual = ForwardDiff.Dual{Float64}(0.0, 1.0)

using KrylovKit
Neig = 8#number of required eigen values
Krylov_arg = (Neig, :LM, KrylovKit.Arnoldi(tol=1e-15, krylovdim=3 + 20, verbosity=0, eager=true));

τmax = τ #maximal timedelay in the mapping
Nstep = 100 # discretization number of the mapping
Timeperiod = Tend # timeperiod of the mapping

#Creating the problem
dp_0_Tfix = dynamic_problemSampled(prob_long, Solver_args, τmax,
    Timeperiod; Historyresolution=Nstep,
    zerofixpont=true, affineinteration=0,
    Krylov_arg=Krylov_arg)

#bv_affine = LinRange(-2.0, 2.05, 201)#TODO: ez kell ha a balolbali szakaszt szeretném követeni
bv_affine = LinRange(2.05, -2.0, 201)
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

Solver_args_T_short = Dict(:alg => MethodOfSteps(RK4()), :verbose => false, :reltol => 1e-6,
    :callback => cb_short)#

Nstep = 50

using ForwardDiff
#TODO: ez nem hiszem, hogy jó megoldás
Base.:convert(::Type{Float64}, x::ForwardDiff.Dual{Float64,Float64,1}) = x.value

using KrylovKit
Neig = 6#number of required eigen values
Krylov_arg = (Neig, :LM, KrylovKit.Arnoldi(tol=1e-30, krylovdim=3 + 35, verbosity=0, eager=true));
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
ustart=saff_c
#bv_affine_H = LinRange(-2.0, -0.951233, 27)
#bv_affine_H = LinRange(-1.0, -4.0, 37)
#bv_affine_H = LinRange(1.2, 0.08, 27)
bv_affine_H = LinRange(b_H_start, -2.0, 150)
bv_affine_H = [LinRange(b_H_start, -1.8, 10)...,LinRange(-1.81, -1.96, 40)...,LinRange(-1.96, -1.98, 100)...]

bv_affine_H = [LinRange(b_H_start, 2, 100)...]
bv_affine_H = [b_H_start:0.01:2.0...]

λ_μ₀_Hopf = Any[similar(bv_affine_H)...]
Amp_H = Any[similar(bv_affine_H)...]
T_period_H = Any[similar(bv_affine_H)...]
#
#plot()
#Threads.@threads
@time for ib in eachindex(bv_affine_H)
    println(ib)
    @show bloc = deepcopy(bv_affine_H[ib])
    @show typeof(bloc)
    mu, saff, sol0 = affine(dp_0_cb, ustart; p=(ζ, δ, bloc, τ, μ))
    #mu, saff, sol0 = affine(dp_Hopf_callback; p=(ζ, δ, bloc, τ, μ))
    λ_μ₀_Hopf[ib] = log.(abs.(mu[1])) / sol0.t[end]
    #Amp[ib] = saff
    T_period_H[ib] = sol0.t[end]
    Amp_H[ib] = maximum(abs.(getindex.(saff, 1)))

    ustart = saff #TODO: enélkül nem áll be szépen (de csak azért nem csinálta, mert affineinteration nem érvényesült)

    #plot!(sol0[1, :], sol0[2, :])
    #aaa=plot!(getindex.(saff, 1), getindex.(saff, 2), marker=:circle, markersize=1, lab="")#
    #display(aaa)


end
#plot!(legend=false)


plot(bv, norm_solperiod)
marcolor=sign.(getindex.(λ_μ₀_Hopf, 1))
scatter!(bv_affine_H, Amp_H,  zcolor=marcolor,color=:bamako,lw=0,markerstrokewidth=0)
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

plot!(bv_affine_H, T_period_H ./ 10.0,lw=3)



## ------------------------------------------------------------------------------

## TODO: csak itt tartok -  tesztelgetés
# Pseudo-archlength method
5 + 5

p=(ζ, δ, b_H_start, τ, μ)

mu_c, saff_c, sol0_c = affine(dp_0_cb, Acrtit_cb * 10.0; p=(ζ, δ, b_H_start, τ, μ));
mu_c_1, saff_c_m1, sol0_c_1 = affine(dp_0_cb, Acrtit_cb * 10.0; p=(ζ, δ, b_H_start-0.01, τ, μ));
scatter(getindex.(saff_c, 1), getindex.(saff_c, 2), lw=2)
scatter!(sol0_c[1, :], sol0_c[2, :])

maximum(abs.(getindex.(saff_c, 1)))
maximum(abs.(getindex.(saff_c_1, 1)))
Δu_Δλ=(saff_c .- saff_c_m1  ) ./ 0.01
#Δu_Δλ=s0 .* 0.0

#mu, saff, sol0 = affine(dp_0_Tfix; p=(ζ, δ, bloc+one_espilon_Dual, τ, μ))
dp=dp_0_cb
Nstep = size(dp.StateSmaplingTime, 1)
s0 = rand(typeof(dp.Problem.u0), Nstep)
s0 .*= 0.0


#initialization of a step:
Δλ=0.01

p=p  .+ (0.0, 0.0, Δλ, 0.0, 0.0)
s0=saff_c+Δu_Δλ*Δλ


## iteration within a step -------

#pDual=p .* 0.0 .+ ForwardDiff.Dual{Float64}(0.0, 1.0)
pDual=p .+ (0.0, 0.0, ForwardDiff.Dual{Float64}(0.0, 1.0), 0.0, 0.0)

v0 = LinMap(dp, s0; p=p)[1]
v0_dual = LinMap(dp, s0; p=pDual)[1]
v0=valuepart.(v0_dual)
dv0dλ=partialpart.(v0_dual)
norm(s0-v0)

TheMapping(s) = partialpart.(LinMap(dp, s * one_espilon_Dual + s0; p=p)[1] - v0)
Nstep = size(dp.StateSmaplingTime, 1)
s_start = rand(typeof(dp.Problem.u0), Nstep)

mus = getindex(schursolve(TheMapping, s_start,dp.Krylov_arg...), [3, 2, 1])


eigval=mus[1]
 eigvec= mus[2]

#a0 = real.(find_fix_pont(s0, v0, mus[1], mus[2],dv0dλ,Δu_Δλ))
x = (v0 - s0)

# println("------------------------------------")
# println(AtA)
 Atx = [dot(eigvec[i], x) for i in 1:size(eigvec, 1)]
 Atdvdlam = [dot(eigvec[i], dv0dλ) for i in 1:size(eigvec, 1)]
 Δu_ΔλAt = [dot(Δu_Δλ,eigvec[i]) for i in 1:size(eigvec, 1)]
 ci =  Atx #TODO: ez ugyan azt adja Schur esetén!!!
 

 error("itt tartok, valahol itt kell előjeleket kereseni talán")
T_jac=vcat(hcat(diagm(eigval .- 1.0),Atdvdlam),
hcat(Δu_ΔλAt',1.0))
@time ci_arch=T_jac\ vcat(Atx,0.0)
 ci=ci_arch[1:end-1]
 Δλ=real.(ci_arch[end])

 ci_mu_original = ci .* (eigval) 
 ci_mu = (Atx .* ((eigval) ./ (eigval .- 1.0)))#TODO: Szabad ezt csinálni, a Schur-nál, nem a sajátértékkel kellenen skálázni... (vagy az pont kiesik valós függvényeknél???)
 #ci_mu-ci_mu_original
 
 #A=transpose(mapreduce(permutedims, vcat, eigvec))
 #fix_v = v0 - A * (((A'A) \ (A' * x)) .* ((eigval) ./ (eigval .- 1.0)))

 fix_v = v0 - mapreduce(x -> x[1] * x[2], +, zip(eigvec, ci_mu))
 s0=real.(fix_v)
 
 @show p=p .+ (0.0, 0.0, Δλ, 0.0, 0.0)
 @show norm(s0-LinMap(dp, s0; p=p)[1])
5+5





































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