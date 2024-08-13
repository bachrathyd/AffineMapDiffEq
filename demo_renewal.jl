using Plots
plotly()


dt=0.001
#history
h(t)=sin.(t*10.0).-sign.(t.+1)+2.0
u0=1.0
u0=h(0.0)

t0=collect(-2.0:dt:0.0)
x0=h.(t0)
x0[end]=u0

#parameters
#a=0.2;#nice
a=1.0#ugly
b=1

x_renewal(t) = a * x(t - 1) + b * x(t - √2).^0.5
x_renewal_callback(t)=x_renewal(t)>4.0 ?  x_renewal(t)/2.0 : x_renewal(t)
x(t)= x0[findfirst(x->x>t, t0)] #interpolation: next closest (ZOH)


T=20
Nstep=ceil(Int,T /   dt)
@time for i in 1:Nstep 
    tnew=t0[end]+dt;
    #push!(x0,x_renewal(tnew))
    push!(x0,x_renewal_callback(tnew))
    push!(t0,tnew)
end
plot(t0,x0)

#---------------------- using neutral DDE solver ----------------
using DifferentialEquations
function renewal(u, h, p, t)
    a,b=p
    du=a * h(p, t - 1.0, Val{1})+ b * h(p, t - √2, Val{1}).^0.5
    du=du>4.0 ? du/2.0 : du
end

p=a,b
h(p, t) = 0.0
h(p, t, deriv::Type{Val{1}}) = t==0.0 ? u0 : h(t)

prob_renewal = DDEProblem(renewal,h, (0.0, T ), p, neutral=true)

@time sol = solve(prob_renewal, MethodOfSteps(Tsit5()),reltol=1e-18)

#plot!(sol) #irrelevant variable
plot!(sol.t,sol(sol.t, Val{1}).u)




error("Ezt itt még nem megy!")

# ---------------- Affine mapping ---------------------------
using DDE_mapping

using KrylovKit
Neig=4#number of required eigen values
Krylov_arg=(Neig,:LM, KrylovKit.Arnoldi(tol=1e-12,krylovdim=8+5,verbosity=3));

τmax=√2 #maximal timedelay in the mapping
Nstep = 100 # discretization number of the mapping
Timeperiod=2.0 # timeperiod of the mapping

Solver_args = Dict(:alg => MethodOfSteps(Tsit5()), :verbose => false, :reltol => 1e-10)#
#Creating the problem
dprenewal = dynamic_problemSampled(prob_renewal, Solver_args, τmax,
Timeperiod; Historyresolution=Nstep,
    zerofixpont=false,    affineinteration=1,
    Krylov_arg=Krylov_arg)


#solvig the problem (mu: Floquet multiplier, saff: discretized foxed point (in [-τmax,0]), sol: the corresponding periodic solution of the fixpoint)
@time mu, saff, sol0 = affine(dprenewal; p=p);


# Comparing the solutions:
plot(sol_period[1, :], sol_period[2, :], marker=:circle,markersize=6,lab="")
plot!(getindex.(saff,1), getindex.(saff,2), marker=:circle,markersize=4,lab="")#
plot!(sol0[1, :], sol0[2, :], marker=:circle,lw = 0,lab="")#marker=:cross,markersize=2)#
##


# Plotting the Floquet multipliers
#in log scale
plot(log.(abs.(mu[1])))
#in complex plane
scatter(mu[1])
plot!(sin.(0:0.01:2pi), cos.(0:0.01:2pi))
