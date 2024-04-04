
includet("src\\DDE_mapping_types.jl")
includet("src\\DDE_mapping_functions.jl")
includet("src\\fucntion_composition.jl")


using Plots



f0 = funcomp(cos, [1.0, 2pi])
function Operator(his::funcomp{Float64,Float64})::funcomp{Float64,Float64}# where Tout  #{T,Tout}

    #funcomp(1.0, (t)->his(t)/t, his.range)
    funcomp(1.0, (t)->his(t+0.25), his.range)
end

#-------------Linearitás teszt-----------
#f0 = funcomp((t)->log(t)+1.0, [1.0, 2pi])
#f1 = funcomp((t)->sin(t),[1.0, 2pi])
f0 = funcomp(cos, [1.0, 2pi])
f1 = funcomp(sin,[1.0, 2pi])
f12=funcomp([f0,f1],[1.0, 2pi])


fout0=Operator(f0)
fout1=Operator(f1)
fout12=Operator(f12)


plot(t,f0.(t)+f1.(t))
#plot!(t,f1.(t))
plot!(t,f12.(t))
plot!(t,fout0.(t) .+ fout1.(t))
#plot!(t,fout1.(t))
plot!(t,fout12.(t))
#-------------Linearitás teszt-----------


# <<<<<<<<<<< Többszörös leképezés teszt >>>>>>>>>>>>>

fout=Operator(f0)

t=fout.range[1]:0.01:fout.range[2]
plot(t,fout.(t))
for _ in 1:20
    fout=Operator(fout)
    plot!(t,fout.(t))
end
plot!()


# <<<<<<<<<<< Többszörös leképezés teszt >>>>>>>>>>>>>

f0 = funcomp(t->1.0+t^2, [1.0, 2pi])
musFun = getindex(schursolve(Operator, f0, 3, :LM,
        KrylovKit.Arnoldi(krylovdim=10, tol=1e-5, verbosity=3)), [3, 2, 1, 4]);


#------------- eigen functions -----------------
NN = 1#1:2:7
plot(t, getindex.([musFun[2][NN](tt) for tt in t], 1), width=2, linestyle=:dash,lab=abs(musFun[1][NN]))
for NN in 2: length(musFun[2])
    plot!(t, getindex.([musFun[2][NN](tt) for tt in t], 1), width=2, linestyle=:dash,lab=abs(musFun[1][NN]))
end
plot!()
#------------- eigen functions -----------------
NN = 1#1:2:7
plot(t, getindex.([musFun[2][NN](tt)*real(musFun[1][NN]) for tt in t], 1), width=2, linestyle=:dash,lab="shape")
#Mapped:
fmapped=Operator(musFun[2][NN])
plot!(t, getindex.([fmapped(tt).*real(musFun[1][NN]) for tt in t], 1), width=2, linestyle=:dash,lab="mapped")


