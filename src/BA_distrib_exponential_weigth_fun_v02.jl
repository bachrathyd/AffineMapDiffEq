#Bartfai Andras> ENOC:turning with acceleration feedback
# citation:
# Single degree of freedom orthogonal cutting model subjected to
# distributed acceleration feedback control
# Andras Bartfai · Zsolt Iklodi · Zoltan Dombovari
5 + 5


#using Revise
using DDE_mapping

using BenchmarkTools
using Plots
plotly()
#gr()
using Profile
using StaticArrays
using DifferentialEquations

# using Memoization

#using MDBM


using Integrals

function diff_cutting_acc_control_neutral(u, h, p, t)
    w,Ωrpm, ωₙ, ζ, μ, k,ρ1,ρ2,ρ3,η,ηnl , τN, τ1,σ1,σ2,σ3,σ4,σ5= p
    
    τc=60/Ωrpm;
    #state variables
    xt = u[1]
    vt = u[2]
    #delayed positon
    xt_τc= h(p, t - τc)[1]
    
    #chip thickness
    hdin=h0+xt_τc-xt
    #cutting force
    Fx = w * (ρ1 * hdin + ρ2 * hdin^2 + ρ3 * hdin^3)
    
    #\\TODO: ez csak egy gyors test a distributed-re, de ettől nagyon lassú lesz, és még le kell tesztelni!
    #distributes
    at_τ_foo(t, p) = h(p, t, Val{1})[2]# derivative of the delayed velocity --> delayed acceleration
    ker_lin(θ, p) =at_τ_foo(t+θ, p)*η(θ)
    ker_cub(θ, p) =at_τ_foo(t+θ, p)^3*ηnl(θ)
    prob_lin = IntegralProblem(ker_lin,  (- τN, - τ1),p)
    prob_cub = IntegralProblem(ker_cub, (- τN, - τ1),p)
    #Int_reltol=1e-7
    #Int_abstol=1e1
    solInt_lin = solve(prob_lin, HCubatureJL(); reltol=Int_reltol)#, abstol=Int_abstol)
    solInt_cub = solve(prob_cub, HCubatureJL(); reltol=Int_reltol)#, abstol=Int_abstol)
    #@show solInt_lin.u
    #@show solInt_cub.u
    Fc= solInt_lin.u+solInt_cub.u


    dxdt = vt
    dvdt = -(2 * ζ * ωₙ * vt + ωₙ^2 * xt + μ * xt^3) + ωₙ^2 / k * (Fc + Fx)
    return SA[dxdt, dvdt]

end


function diff_cutting_acc_control_neutral_perturbed(u, h, p, t)
    w,Ωrpm, ωₙ, ζ, μ, k,ρ1,ρ2,ρ3,η,ηnl , τN, τ1,σ1,σ2,σ3,σ4,σ5 = p
    
    τc=60/Ωrpm;
    #state variables
    ut = u[1]
    dudt = u[2]
    #delayed positon
    ut_τc= h(p, t - τc)[1]
    

    #\\TODO: ez csak egy gyors test a distributed-re, de ettől nagyon lassú lesz, és még le kell tesztelni!
    #distributes
    at_τ_foo(t, p) = h(p, t, Val{1})[2]# derivative of the delayed velocity --> delayed acceleration
    ker_lin(θ, p) =at_τ_foo(t+θ, p)*η(θ)
    ker_cub(θ, p) =at_τ_foo(t+θ, p)^3*ηnl(θ)
    prob_lin = IntegralProblem(ker_lin,  (- τN, - τ1),p)
    prob_cub = IntegralProblem(ker_cub, (- τN, - τ1),p)
    #Int_reltol=1e-7
    #Int_abstol=1e1
    solInt_lin = solve(prob_lin, HCubatureJL(); reltol=Int_reltol)#, abstol=Int_abstol)
    solInt_cub = solve(prob_cub, HCubatureJL(); reltol=Int_reltol)#, abstol=Int_abstol)
    #@show solInt_lin.u
    #@show solInt_cub.u
    Fc= solInt_lin.u+solInt_cub.u

    dudt_out = dudt
    dduddt_out = -(2 * ζ * ωₙ * dudt + ωₙ^2 * ut + μ * ut^3) + ωₙ^2 / k *Fc-
    ωₙ^2 / k *w*( σ1*(ut−ut_τc) +
    σ2*(ut^2 +ut_τc^2) +
    σ3*(ut^3 -ut_τc^3)+
    σ4*(ut*ut_τc)+
    σ5*(ut*ut^2-ut^2*ut_τc) )

    #println(SA[dudt_out, dduddt_out])
    #println(t)
    return SA[dudt_out, dduddt_out]

end

Base.:+(a::SVector, b::Bool) = a .+ b


ωₙ=150*2*pi
ζ=0.01
k=20e6
h0=0.1e-3
ρ1 = 6109.6*1e6;#6109.6 N/mm2
ρ2 = −54141.6*1e9;# −54141.6 N/mm3
ρ3 =203769 *1e12;# 203769 N/mm4
σ1 = ρ1 +2*h0 *ρ2 +3* h0^2*ρ3
σ2 = ρ2 +3*h0 *ρ3
σ3 = ρ3
σ4 =  −2*ρ2 −6*h0* ρ3;
σ5 = 3*ρ3


κ=20#kg
κnl=-0.01
τN=100e-3
τ1 =10e-3
∆T = τN −τ1
#η(θ)::Float64=κ/∆T#TODO: wring parametrization of the input
#ηnl(θ)::Float64=κnl/∆T#TODO: wring parametrization of the input
η(θ)::Float64=κ/∆T*exp((θ+τ1))
ηnl(θ)::Float64=κnl/∆T*exp((θ+τ1))#TODO: wring parametrization of the input
  
μ=0;#knl/m

w=1e-3
Ωrpm=7000#8000 rpm
  

pars= (w,Ωrpm, ωₙ, ζ, μ, k,ρ1,ρ2,ρ3,η,ηnl , τN, τ1,σ1,σ2,σ3,σ4,σ5)

Int_reltol=1e-4
Int_abstol=1e1

h(p, t::Float64) = SA[1.0e-5, 0.0]
h(p, t::Float64, deriv::Type{Val{1}}) = SA[0.0, 1.0e-5]
u0 = SA[-1.0e-5, 0.0]

@time diff_cutting_acc_control_neutral_perturbed(u0, h, pars, 0.1);
# @code_warntype diff_cutting_acc_control_neutral_perturbed(u0, h, pars, 0.1)

#test:
Solver_args = Dict(:alg => MethodOfSteps(RK4()), :verbose => false, :reltol => 1e-3)#

Tsim=0.8
prob_DDE = DDEProblem(diff_cutting_acc_control_neutral_perturbed, u0, h, (0.0, Tsim), pars; neutral=true)
@time sol = solve(prob_DDE; Solver_args...);

plot(sol)
# plot(sol(sol.t[end] .- (0.0:0.01:τ*1.0)))


## ---------------------- Spectrum test brute-force--------------------
#scatter()
gr()
using KrylovKit
N_for=158;#2 min
Ωrpm=7000.0

Int_reltol=1e-4
Int_abstol=1e1

Solver_args = Dict(:alg => MethodOfSteps(RK4()), :verbose => false, :reltol => 1e-5)#
#Solver_args = Dict(:alg => MethodOfSteps(BS3()), :verbose => false, :reltol => 1e-3)#
τN = 20e-3#LinRange(11e-3,30e-3,N_for+1)# initial grid in x direction
τNv = LinRange(11e-3,100e-3,N_for)# initial grid in x direction
τNv = LinRange(11e-3^0.5,150e-3^0.5,N_for).^2# initial grid in x direction


N_for=200#240;
τNv =LinRange(11e-3^0.5,100e-3^0.5,N_for).^2# initial grid in x direction
TN=5
#TNv=LinRange(1,100,N_for)
NTper=20;
#NTperv=(1:N_for)
kiter=N_for
@time for kiter in 1:N_for # initial grid in x direction
    τN=τNv[kiter]
    #TN=TNv[kiter]
    #NTper=NTperv[kiter]
Neig=12#number of required eigen values
Krylov_arg=(Neig,:LM, KrylovKit.Arnoldi(tol=1e-8,krylovdim=Neig+80+kiter÷2,verbosity=1,maxiter = 30));

τmax =maximum([60/Ωrpm,τN,τ1])
Tnatural=2pi/ωₙ 
dtsampling=Tnatural /TN

@show Nstep = Int(τmax ÷ Tnatural*100)

#N=4
#Tperiod=τmax/N
Tperiod=Tnatural/NTper
#Creating the problem
dpdp = dynamic_problemSampled(prob_DDE, Solver_args, τmax,
Tperiod; Historyresolution=Nstep,
    zerofixpont=true,    affineinteration=0,
    Krylov_arg=Krylov_arg)

# fix point by affine map
ploc=deepcopy((w,Ωrpm, ωₙ, ζ, μ, k,ρ1,ρ2,ρ3,η,ηnl , τN, τ1,σ1,σ2,σ3,σ4,σ5))
@time mu, saff, sol0 = affine(dpdp,p=ploc);

#using Profile
#@profview affine(dpdp,p=ploc);

#plot(log.(abs.(mu[1])))
#scatter((mu[1]))
#plot!(sin.(0:0.01:2pi), cos.(0:0.01:2pi))


lams=log.(mu[1])/Tperiod
#aaa=scatter!(τN .* ones(size(lams,1)),real.(lams),legend = false)
#aaa=scatter!(ylim=(-100,100),xlabel="τN",ylabel="Re(λ)")

aaa=scatter(lams,legend = false,xlabel="Re(λ)",ylabel="Im(λ)"
,xlim=(-100,50),ylim=(-4000,4000),title="τN: $τN ")

#prob_DDE = DDEProblem(diff_cutting_acc_control_neutral_perturbed, u0, h, (0.0, Tsim), ploc; neutral=true)
@time sol = solve(remake(prob_DDE,p=ploc,tspan=(0,0.5)); Solver_args...);
bbb=plot(sol,xlabel="t",ylabel="u,v",title="time sim")

plot(dpdp.StateSmaplingTime,getindex.(mu[2][1],1))
ccc=plot!(dpdp.StateSmaplingTime,getindex.(mu[2][1],2),title="Schur mode 1")
plot(dpdp.StateSmaplingTime,getindex.(mu[2][3],1))
ddd=plot!(dpdp.StateSmaplingTime,getindex.(mu[2][3],2),title="Schur mode 3")
figout=plot(aaa,bbb,ccc,ddd,size=(1600,1000))
display(figout)

savefig("test$kiter.pdf")
end
#scatter!(ylim=(-15.5,-14.5))
#TODO: a sajátvektor tesztelése szuimulációval, hogy tényleg elszáll-e
v0,sol_Ai=DDE_mapping.LinMap(dpdp, mu[2][1]; p=ploc)

plot(sol_Ai)
plot!(dpdp.StateSmaplingTime,getindex.(mu[2][1],1),xlim=(dpdp.StateSmaplingTime[1], dpdp.Tperiod))
plot!(dpdp.StateSmaplingTime,getindex.(mu[2][1],2),title="Schur mode 1")
plot!(dpdp.StateSmaplingTime .+ dpdp.Tperiod,getindex.(v0,1))
plot!(dpdp.StateSmaplingTime.+ dpdp.Tperiod,getindex.(v0,2),title="Schur mode 1",xlim=(dpdp.StateSmaplingTime[1], dpdp.Tperiod))

## --------------------

println("----------Start brute-force---------------")

Neig=2#number of required eigen values
Krylov_arg=(Neig,:LM, KrylovKit.Arnoldi(tol=1e-4,krylovdim=Neig+8,verbosity=0));

N=100 #12 min
N=60;#2 min
Ωrpmv = LinRange(2000,9000,N+1)# initial grid in x direction
wv = LinRange(-0.1e-3,10e-3,N)# initial grid in y direction
Spek = zeros(size(wv, 1), size(Ωrpmv, 1))
#Threads.@threads

@time Threads.@threads  for j in 1:size(Ωrpmv, 1)
    println(j)
    Ωrpm = Ωrpmv[j]

    Nstep = 40
    N=20
    τmax =maximum([60/Ωrpm,τN])
    
    T=τmax/N
    #Creating the problem
    dploc = dynamic_problemSampled(prob_DDE, Solver_args, τmax,
    T; Historyresolution=Nstep,
        zerofixpont=true,    affineinteration=0,
        Krylov_arg=Krylov_arg)
        Threads.@threads for i in 1:size(wv, 1)
        #println([i,j])
        w = wv[i]
        mu, saff, sol0 = affine(dploc; p=(w,Ωrpm, ωₙ, ζ, μ, k,ρ1,ρ2,ρ3,η,ηnl , τN, τ1,σ1,σ2,σ3,σ4,σ5))
        muMAX = abs(mu[1][1].^N)
        #muMAX = spectralradius(dploc; p=(ωn, τf, τr, w, μ))
        Spek[i, j] = muMAX
    end
end

Spek_sat = deepcopy(Spek);
Spek_sat[Spek_sat.>1.0] .= 1.0;
heatmap(Ωrpmv,wv, (Spek_sat),xlabel="w", ylabel="Ωrpm")




## --------------------

println("----------Start brute-force---------------")

Neig=2#number of required eigen values
Krylov_arg=(Neig,:LM, KrylovKit.Arnoldi(tol=1e-4,krylovdim=Neig+8,verbosity=0));

N=10;#2 min
τNv = LinRange(0.01e-3,200e-3,N+1)# initial grid in x direction
wv = LinRange(-0.1e-3,10e-3,N)# initial grid in y direction
Spek = zeros(size(wv, 1), size(τNv, 1))
#Threads.@threads

@time Threads.@threads  for j in 1:size(τNv, 1)
    println(j)
    Ωrpm = 7000
    τN=τNv[j]

    Nstep = 40
    N=20
    τmax =maximum([60/Ωrpm,τN])
    
    T=τmax/N
    #Creating the problem
    dploc = dynamic_problemSampled(prob_DDE, Solver_args, τmax,
    T; Historyresolution=Nstep,
        zerofixpont=true,    affineinteration=0,
        Krylov_arg=Krylov_arg)
        Threads.@threads for i in 1:size(wv, 1)
        #println([i,j])
        w = wv[i]
        mu, saff, sol0 = affine(dploc; p=(w,Ωrpm, ωₙ, ζ, μ, k,ρ1,ρ2,ρ3,η,ηnl , τN, τ1,σ1,σ2,σ3,σ4,σ5))
        muMAX = abs(mu[1][1].^N)
        #muMAX = spectralradius(dploc; p=(ωn, τf, τr, w, μ))
        Spek[i, j] = muMAX
    end
end

Spek_sat = deepcopy(Spek);
Spek_sat[Spek_sat.>1.0] .= 1.0;
heatmap(τNv,wv, (Spek_sat),xlabel="w", ylabel="τN")



##----------------------  stability map --------------------
#using MDBM
#
##ax1 = Axis(-0.009:0.01:0.0005, "τf") # initial grid in x direction
##ax2 = Axis(-0.5:0.25:0.5, "K") # initial grid in y direction
#ax1 = Axis(LinRange(-0.009,0.02,15), "τf") # initial grid in x direction
#ax2 = Axis(LinRange(-1.0,1.001651,6), "K") # initial grid in y direction
#function fooMathieu(τf, K)
#    τmax = maximum([τr,τr + τf]) * 1.2
#    Nstep = 200
#    dploc = dynamic_problemSampled(pro_BD_Enoc, Solver_args, τmax,
#T; Historyresolution=Nstep,
#    zerofixpont=true,    affineinteration=0,
#    Krylov_arg=Krylov_arg)
#    w(t) = K*exp(A*t);
#    #w(t) = K
#    #println((τf, K))
#    mu, saff = affine(dploc; p=(ωn, τf, τr, w, μ))
#    
#    ABSmuMax = abs(mu[1][1]) ;#abs(mu[1][3])
#    return ABSmuMax - 1.0
#    
#    #MuMin1_prod = prod(abs.(mu[1][1:4]) .-1 ) ;#abs(mu[1][3])
#    #return MuMin1_prod#ABSmuMax - 1.0
#end
#
#mymdbm = MDBM_Problem(fooMathieu, [ax1, ax2])
#iteration = 3#number of refinements (resolution doubling)
#@time MDBM.solve!(mymdbm, iteration)
##points where the function foo was evaluated
#x_eval, y_eval = getevaluatedpoints(mymdbm)
##interpolated points of the solution (approximately where foo(x,y) == 0 and c(x,y)>0)
#x_sol, y_sol = getinterpolatedsolution(mymdbm)
##scatter(x_eval,y_eval,markersize=1)
##scatter!(x_sol,y_sol,markersize=3)
#scatter!(x_sol, y_sol, color=:blue, markersize=3)
#

#---------------------- BA-D-curve solution --------------------
using MDBM
using Plots


Nstart=8
ax1 = Axis(LinRange(3000.0,9000.0,Nstart), "Ωrpm") # initial grid in y direction
ax2 = Axis(LinRange(-1e-3,8-3,Nstart), "w") # initial grid in x direction
ax3 = Axis(LinRange(-10.0,1500.0,Nstart), "ωc") # initial grid in y direction
function foo_BA_Dcurve(Ωrpm_loc::Float64, w_loc::Float64,ω_loc::Float64)::SVector{2, Float64}


    ωₙ=150*2*pi
    ζ=0.01
    k=20e6
    h0=0.1e-3
    ρ1 = 6109.6*1e6;#6109.6 N/mm2
    ρ2 = −54141.6*1e9;# −54141.6 N/mm3
    ρ3 =203769 *1e12;# 203769 N/mm4
    σ1 = ρ1 +2*h0 *ρ2 +3* h0^2*ρ3
  
    
    κ=20#kg
    τN=10e-3
    τ1 =20e-3
    ∆T = τN −τ1
    η(θ::Float64)::Float64=κ/∆T#TODO: wring parametrization of the input
      

    τc=60/Ωrpm_loc;
    λ=1.0im*ω_loc

    #at_τ_foo(t, p) = h(p, t, Val{1})[2]# derivative of the delayed velocity --> delayed acceleration
    ker_lin(θ::Float64, _)::ComplexF64 =λ^2*exp(λ*θ)*η(θ)
    prob_lin = IntegralProblem(ker_lin,  (- τN, - τ1),[])
    solInt_lin = solve(prob_lin, HCubatureJL(); reltol=1e-4)
    INTsol::ComplexF64=solInt_lin.u

    D=λ^2+2*ζ*ωₙ*λ+  ωₙ^2-ωₙ^2/k*INTsol+ωₙ^2/k*w_loc*σ1*(1-exp(-λ*τc))
    return SA[real(D), imag(D)]::SVector{2, Float64}
end
@time foo_BA_Dcurve(1000.0,5e-3,300.0)
#@benchmark foo_BA_Dcurve(1000.0,5e-3,300.0)
#@code_warntype foo_BA_Dcurve(1000.0,5e-3,300.0)
Dcurve_mdbm = MDBM_Problem(foo_BA_Dcurve, [ax1, ax2, ax3])
iteration = 1#number of refinements (resolution doubling)

#for _ in 1:2
@time MDBM.solve!(Dcurve_mdbm, iteration);
Dcurve_x_sol, Dcurve_y_sol, Dcurve_z_sol = getinterpolatedsolution(Dcurve_mdbm);
@show scatter(Dcurve_x_sol, Dcurve_y_sol, zcolor=Dcurve_z_sol, markersize=2)
@show scatter(Dcurve_x_sol, Dcurve_z_sol, zcolor=Dcurve_z_sol, markersize=2)
#end
