5+5


using Revise
using DDE_mapping

using BenchmarkTools
using Plots
plotly()
#using PlotlyBase
#using PlotlyKaleido
using Profile
using StaticArrays
using DifferentialEquations
function f_now(p, t)
    ζ, δ, ϵ, b, τ , T  = p
   SA[-(δ+ϵ*     cos(t)), -2*ζ]
   #SA[-(δ+ϵ*sign(cos(t))), -2*ζ]#Meissner  
end
function f_past(p, t)
    ζ, δ, ϵ, b, τ , T  = p
    SA[b, 0]
end

function DelayMathieu(u, h, p, t)
    ζ, δ, ϵ, b, τ , T = p
    dx = u[2]
    ddx = f_now(p, t)' * u + b * h(p, t - τ)[1] 
    SA[dx, ddx]
end
Base.:+(a::SVector, b::Bool) = a .+ b


##<<<<<<<<<<< Lin Map based
NF=Float64
ζ = NF(0.02)          # damping coefficient
δ = NF(1.5)#0.2          
ϵ = NF(0.15)#4#5#8;#5      
τ = NF(2pi)          # Time delay
b = NF(0.5)
T= NF(2pi)
p = ζ, δ, ϵ, b, τ , T

u0 = SA{NF}[1.0, 0.0]
h(p, t) = SA{NF}[0.0; 0.0]
probMathieu = DDEProblem(DelayMathieu, u0, h, (NF(0.0), NF(T * 1)), p; constant_lags=[τ])


b = 0.2
p = ζ, δ, ϵ, b, τ, T

Nstepv = floor.(Int, exp.(LinRange(log(3), log(1e5), 20))) #21 # Ez jól néz ki
Tmean = zeros(Float64, size(Nstepv, 1))
Tstd = zeros(Float64, size(Nstepv, 1))
mumaxVec_error = zeros(Float64, size(Nstepv, 1))
BenchmarkTools.DEFAULT_PARAMETERS.samples = 2#1000
BenchmarkTools.DEFAULT_PARAMETERS.seconds = 0.1#5.0
BenchmarkTools.DEFAULT_PARAMETERS.time_tolerance = 0.50#0.1

#the_alg=MethodOfSteps(ABM43())
#the_alg = MethodOfSteps(BS3())
the_alg=MethodOfSteps(RK4())
#the_alg=MethodOfSteps(Tsit5())
τmax = 2pi + 0.1
#dpdp_Time = dynamic_problemSampled(probMathieu, the_alg, τmax, T; Historyresolution=maximum(Nstepv) * 100, eigN=1, zerofixpont=true)
dpdp_Time = dynamic_problemSampled(probMathieu, the_alg, τmax, T; 
Historyresolution=maximum(Nstepv) * 10, eigN=1, zerofixpont=false)
@time mumaxFINE = spectralradius(dpdp_Time; p=p)
for kk in 1:size(Nstepv, 1)
    runpercent = kk / size(Nstepv, 1) * 100
    #Nstep = 100
    Nstep = Nstepv[kk]
    

    dpdp_Time = dynamic_problemSampled(probMathieu, the_alg, τmax, T; Historyresolution=Nstep, eigN=1, zerofixpont=true)

    t = @benchmark spectralradius($dpdp_Time; p=$p)
    muMaxLoc = spectralradius(dpdp_Time; p=p)
    locerror=abs(muMaxLoc - mumaxFINE)
    mumaxVec_error[kk] = locerror
    Tmean[kk] = BenchmarkTools.median(t).time / 1e9
    Tstd[kk] = BenchmarkTools.std(t).time / 1e9

    
    locerror=mumaxVec_error[kk] 
    print("  - run [%]: $runpercent, resolution [1]: $Nstep ,  muerror: $locerror  ")
    print(t)

    println(" <<<<<<< ")
end

#Xax2plot = Nstepv
#XLABEL="resolution [1]"
#Xax2plot=mumaxVec_error
#XLABEL="Mu error"

Xax2plot= Tmean
XLABEL="T-mean [s]"


Yax2plot=mumaxVec_error
YLABEL="Mu error"
#Yax2plot = Nstepv
#YLABEL="resolution [1]"
#Yax2plot = Tmean
#YLABEL="T-mean [s]"

#plot(Xax2plot,Tmean,grid=true,yerror=Tstd)
plot(Xax2plot, Yax2plot, grid=false,  fillalpha=0.5,
    label=false, yaxis=:log10, xaxis=:log10, color=:green)

for kC in (10.0 .^ (-20:20))
    list=[1e-20,1e20]
    plot!(list, kC .* list .^(-0.3333), label=false, yaxis=:log10, xaxis=:log10)
    plot!(list, kC .* list .^(-3.0), label=false, yaxis=:log10, xaxis=:log10)
    plot!(list,kC .* list .^1 , label=false,linestyle=:dash)
    #plot!(list,kC .* list .^2 , label=false,linestyle=:dash)
end

p=plot!(Xax2plot, Yax2plot, color=:green, marker=(:circle, 3, 1.0), label="T-smooth",
    xlabel=XLABEL, ylabel=YLABEL, yaxis=:log10, xaxis=:log10,
    linewidth=3,
    xticks=10.0 .^ (-20:1.0:20), yticks=10.0 .^ (-20:1.0:20),
    xlims=(-maximum(-Xax2plot) / 3, maximum(Xax2plot) * 3),#xlims=(1e-10,1e1),#
    ylims=(-maximum(-Yax2plot) / 3, maximum(Yax2plot) * 3),
    grid=true,
    gridlinewidth=2)

#p=plot!( xlims=(1e-0,1e5),ylims=(1e-6,1e0))
#p=plot!( xlims=(1e-15,1e2),ylims=(1e-0,1e4))
  #  display(plt), gui()
   default(show = true)
   default(show = false)
#savefig(p, "smooth_non_smooth_Mathieu__fix_Timestep.png")
#savefig("smooth_non_smooth_Mathieu__fix_Timestep.png") 