#OpeartorMap tests
5+5

using LinearAlgebra
 
using StaticArrays
using DifferentialEquations
 
using Profile
using BenchmarkTools
 
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 5
BenchmarkTools.DEFAULT_PARAMETERS.samples = 100
 
using Plots
plotly()
 
### -------------- Mathieu solution to have a function like stuff -----------
function DelayMathieu(u, h, p, t)
    ζ, δ, ϵ, b, τ, T = p  # Parameters
    F = 0.1 * (cos(2pi * t / T) .^ 10)#External forcing
    # Components of the delayed differential equation
    dx = u[2]
    ddx = -(δ + ϵ * cos(2pi * t / T)) * u[1] - 2 * ζ * u[2] + b * h(p, t - τ)[1] + F
    # Update the derivative vector
    SA[dx, ddx]
    #MVector(dx, ddx)
end
 
## parameters
ζ = 0.02          # damping coefficient
δ = 1.5#0.2          # nat. freq
ϵ = 0.15#4#5#8;#5         # cut.coeff
τ = 2pi          # Time delay
b = 0.5
T = 2pi
p = ζ, δ, ϵ, b, τ, T
#p = (ζ, ωn, k, τ,10.0)
 
 
 
# test simulation ---------------
#history function
h(p, t) = SA[10.0; 0.0]
#h(p, t) = MVector(10.0, 0.0)
##initial condition
#u0 = SA[1.0, 0.0]
#probMathieu = DDEProblem(DelayMathieu, u0, h, (0.0, T * 1000.0), p; constant_lags=[τ])
 
Tmapping = 2π
probMathieu = DDEProblem(DelayMathieu, h, (0.0, Tmapping), p; constant_lags=[τ])
#Parameters for the solver as a Dict (it is necessary to collect it for later use)
Solver_args = Dict(:alg => MethodOfSteps(RK4()), :verbose => false, :reltol => 1e-4)#
 
##TODO: maybe better with integrator "Integrator Interface" integrator = init(prob, alg; kwargs...) - it contains the necessary info for the solver
#integrator_sol = init(probMathieu; Solver_args...)
#solve!(integrator_sol)
#plot(integrator_sol.sol)
#integrator_sol.opts
#integrator_sol.alg
sol = solve(probMathieu; Solver_args...)
plot(sol)

#@benchmark solve(probMathieu; Solver_args...)
 
## 
 

# ---------------- Affine mapping ---------------------------
using Revise
using DDE_mapping
using KrylovKit #TODO: Now it is @6.1 version, ugrade by: pin DataFrames leading to errors

Base.:+(a::SVector, b::Bool) = a .+ b
Base.:+(a::SVector, b::Float64) = a .+ b #TODO: where to put this?

Base.:+(a::MVector, b::Bool) = a .+ b
Base.:+(a::MVector, b::Float64) = a .+ b #TODO: where to put this?


Neig=4#number of required eigen values
Krylov_arg=(Neig,:LM, KrylovKit.Arnoldi(tol=1e-8,krylovdim=Neig+6,verbosity=0));

τmax=τ #maximal timedelay in the mapping
Nstep = 100 # discretization number of the mapping
Timeperiod=Tmapping # timeperiod of the mapping

#Creating the problem
dpMathieu = dynamic_problemSampled(probMathieu, Solver_args, τmax,
Timeperiod; Historyresolution=Nstep,
    zerofixpont=false,    affineinteration=1,
    Krylov_arg=Krylov_arg)

@time mu, saff, sol0 = affine(dpMathieu; p=p);
#  @benchmark affine(dpMathieu; p=p)
   
  # Comparing the solutions:
  #plot(sol_period[1, :], sol_period[2, :], marker=:circle,markersize=6,lab="")
  plot(getindex.(saff,1), getindex.(saff,2), marker=:circle,markersize=4,lab="")#
  plot!(sol0[1, :], sol0[2, :], marker=:circle,lw = 0,lab="")#marker=:cross,markersize=2)#
  ##
  
  
  # Plotting the Floquet multipliers
  #in log scale
  plot(log.(abs.(mu[1])))
  #in complex plane
  scatter(mu[1])
  plot!(sin.(0:0.01:2pi), cos.(0:0.01:2pi))
  
## ----------testing the operator - like behavioure
 
 

using FunctionWrappers
import FunctionWrappers: FunctionWrapper
 
using Revise
include("fucntion_composition.jl")
#include("sol_operator_functions.jl")
include("ReangedFun.jl")

mutable struct CallbackS2_F64 <: Function
    f::FunctionWrapper{SVector{2, Float64},Tuple{Float64}}
end
(cb::CallbackS2_F64)(v) = cb.f(v)

fXX = CallbackS2_F64((t) -> sol.prob.h(p,t))
typeof(fXX)
#@code_warntype fXX(3.1415)
#fXX(3.1415)


FX_FC=funcomp(fXX,[-τ,0.0])
#FX_FC=funcomp{Float64,SVector{2, Float64}}((t) -> sol.prob.h(p,t),[-τ,0.0])
FX_FC_comp=funcomp((t::Float64) -> sol.prob.h(p,t)::SVector{2, Float64},[-τ,0.0])
FX_FC_Range=RangeFun{Float64,SVector{2, Float64}}((t) -> sol.prob.h(p,t),[-τ,0.0])


# ~~~~~~~~~~~~~~~~~~~~~ RangeFunction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

foo=FunctionWrapper{SVector{2, Float64},Tuple{Float64}}( (t) -> sol.prob.h(p,t) )
aaa=RangeFun{Float64,SVector{2, Float64}}(foo,MVector(-τ,0.0))

aaa(1.2)
sol.prob.h(p,1.2)
ccc=aaa-aaa
typeof(aaa)
aaa(1.0)
bbb= similar(aaa)
 #@code_warntype
 aaa(1.0)
 empty!(bbb)
 typeof(bbb)
 bbb(1.1)

 LinearAlgebra.mul!(bbb,aaa,5.0)
 bbb(1.1)
 aaa(1.1)
 LinearAlgebra.rmul!(bbb,0.1)
 LinearAlgebra.axpy!(10.0,aaa,bbb)
 LinearAlgebra.axpby!(10.0,aaa,0.1,bbb)
 5+5
 bbb(1.1)
 #scale!!(bbb,0.1)
 bbb(1.1)
# ~~~~~~~~~~~~~~~~~~~~~ RangeFunction ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#function Operator_Sol(sol_in)#::T where T
#    Tmap=sol_in.t[end]
#   
#    fX_sol = CallbackS2_F64(sol_in)
#    fX_h= CallbackS2_F64((tloc)->sol_in.prob.h(sol_in.prob.p,tloc))
#    u0=sol_in.prob.u0
#    function h_hist(p,tin)
#        t=(tin-Tmap)
#        if t < 0.0
#           return fX_h(t)::typeof(sol.prob.u0)
#        elseif t == 0.0
#            u0::typeof(sol.prob.u0)
#        else
#            fX_sol(t)::typeof(sol.prob.u0)
#        end
#    end
#     sol_out = solve(remake(sol_in.prob; u0=sol_in(Tmap), tspan=(0.0, Tmap), h=h_hist);alg = MethodOfSteps(RK4()), verbose = false, reltol = 1e-4)
#   return sol_out
#end
 
 

function Operator_FunX(foo::CallbackS2_F64,dp::dynamic_problemSampled)
    sol_out = solve(remake(dp.Problem; u0=foo(0), tspan=(0.0, dp.Tperiod), h=(p,tin)->foo(tin)); dp.alg...)
    fX_sol_out = CallbackS2_F64((t_find)->DDE_mapping.getvalues(sol_out,t_find+ dp.Tperiod))
    return fX_sol_out::CallbackS2_F64 ,sol_out
end

function Operator_FunX(foo::funcomp,dp::dynamic_problemSampled)
    #println("-----------zzzzzzzzzzzzzzzzzzzzzzz-------")
    #println(foo)
   #@time
    sol_out = solve(remake(dp.Problem; u0=foo(0.0), tspan=(0.0, dp.Tperiod), h=(p,tin)->foo(tin)); dp.alg...)
    fX_sol_out = funcomp((t_find)->DDE_mapping.getvalues(sol_out,t_find+dp.Tperiod),[-dp.maxdelay,0.0])
    #println("-----------<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>-------")
    return fX_sol_out::funcomp,sol_out
end

function Operator_FunX(foo::RangeFun{Tin,Tout},dp::dynamic_problemSampled) where {Tin,Tout}
    sol_out = solve(remake(dp.Problem; u0=foo(0.0), tspan=(0.0, dp.Tperiod), h=(p,tin)->foo(tin)); dp.alg...)
    fX_sol_out = RangeFun{Tin,Tout}((t_find)->DDE_mapping.getvalues(sol_out,t_find+dp.Tperiod),[-dp.maxdelay,0.0])
    return fX_sol_out::RangeFun{Tin,Tout},sol_out
end



fXX(-1.0)
FX_FC(-1.0)
f0=FX_FC_Range
f0(1.0)
f0=FX_FC_comp
@time f1=Operator_FunX(f0,dpMathieu)[1];
@time f2=Operator_FunX(f1,dpMathieu)[1];
@time f2=Operator_FunX(f1,dpMathieu)[1];
#@benchmark Operator_FunX(f2,dpMathieu)




#  long single simulation 
Ns=100
probMathieu_long = DDEProblem(DelayMathieu, h, (0.0, Ns*Tmapping), p; constant_lags=[τ])
@time sol_long = solve(probMathieu_long; Solver_args...);
#@benchmark sol_long = solve(probMathieu_long; Solver_args...)
plot(sol_long)

#long simulation by applying the Operator repeteadly
t_past=LinRange(-τ,0,200)
s=FX_FC
s=FX_FC_Range
kk=0
v,soliter=Operator_FunX(s,dpMathieu);
kk +=1
#plot(sol)
plot!(t_past .+ (kk-1)*dpMathieu.Tperiod,getindex.(s.(t_past),1))#,xlim=[-τ,kk*dpMathieu.Tperiod])
plot!(soliter.t.+(kk-1)* dpMathieu.Tperiod,getindex.(soliter.u,1))
plot!(t_past .+(kk)* dpMathieu.Tperiod,getindex.(v.(t_past),1))#,xlim=[-τ,kk*dpMathieu.Tperiod])
s=v

#TODO: ez kb 3 - 5 ször lassab (plottolás nélkül), de ezt talán annyira nem gond
@time for _ in 1:Ns

    v,soliter=Operator_FunX(s,dpMathieu);
    kk +=1
   plot!(t_past .+ (kk-1)*dpMathieu.Tperiod,getindex.(s.(t_past),1))#,xlim=[-τ,kk*dpMathieu.Tperiod])
   plot!(soliter.t.+(kk-1)* dpMathieu.Tperiod,getindex.(soliter.u,1))
   plot!(t_past .+(kk)* dpMathieu.Tperiod,getindex.(v.(t_past),1))#,xlim=[-τ,kk*dpMathieu.Tperiod])
    
    s=v
end
plot!()








## -------------- KrylovKit Operator testing -----------------------
FX_FC_comp
FX_FC_Range
eltype(FX_FC_comp)
similar(FX_FC_comp,Int)

scalartype(FX_FC_comp)
scalartype(typeof(FX_FC_comp))

scalartype(FX_FC_comp)
norm(FX_FC_comp)
dot(FX_FC_comp,FX_FC_comp)
scalartype(FX_FC_Range)
norm(FX_FC_Range)
dot(FX_FC_Range,FX_FC_Range)






#scale!!(FX_FC_Range,0.1)
#scale!!(FX_FC_comp,0.1)
s0=FX_FC_Range;
s0=FX_FC_comp;
v0 = Operator_FunX(s0,dpMathieu)[1];
#println(norm(s0-v0))
s_start =   Operator_FunX(v0,dpMathieu)[1];
 
s0(1.0)
#scale!!(s0,0.1);
schursolve((xx)->Operator_FunX(xx+s0,dpMathieu)[1]-v0, s_start,Krylov_arg...)

@time musOP = getindex(schursolve((xx)->Operator_FunX(xx+s0,dpMathieu)[1]-v0, s_start,Krylov_arg...), [3, 2, 1]);
#@show @benchmark getindex(schursolve((xx)->Operator_FunX(xx+s0,dpMathieu)[1]-v0, s_start,Krylov_arg...), [3, 2, 1])



# Plotting the Floquet multipliers
#in log scale
plot(log.(abs.(mu[1])))
plot!(log.(abs.(musOP[1])))
#in complex plane
scatter(mu[1])
scatter!(musOP[1])
plot!(sin.(0:0.01:2pi), cos.(0:0.01:2pi))


#musOP = eigsolve((xx)->Operator_FunX(xx+s0,dpMathieu)[1]-v0, s_start, 8, :LM)
#
#
#t_past=LinRange(-τ,0,200)
#A1=reduce!(musOP[2][1])
#A1map,solA1=Operator_FunX(A1,dpMathieu)
#plot(t_past,getindex.(A1.(t_past),1))
#plot!(t_past,getindex.(A1map.(t_past),1))
#
#plot!(t_past,real.(getindex.(musOP[1][1]*A1.(t_past),1)))
#musOP[1][1]
#
#musOP[1][1]*A1.(t_past)-A1map.(t_past)
#




