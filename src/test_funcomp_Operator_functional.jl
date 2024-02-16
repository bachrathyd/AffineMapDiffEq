
includet("src\\DDE_mapping_types.jl")
includet("src\\DDE_mapping_functions.jl")
includet("src\\fucntion_composition.jl")


using Plots


#-------------Linearitás teszt-----------
#f0 = funcomp((t)->log(t)+1.0, [1.0, 2pi])
#f1 = funcomp((t)->sin(t),[1.0, 2pi])
f0 = funcomp(cos, [1.0, 2pi])
f1 = funcomp(sin,[1.0, 2pi])
f12=funcomp([f0,f1])


fout0=Operator(f0)
fout1=Operator(f1)
fout12=Operator(f12)


plot(t,f0.(t)+f1.(t))
#plot!(t,f1.(t))
plot!(t,f12.(t))
plot!(t,fout0.(t) .+ fout1.(t))
#plot!(t,fout1.(t))
plot!(t,fout12.(t))


function Operator(his::funcomp{Float64,Float64})::funcomp{Float64,Float64}# where Tout  #{T,Tout}
 #   function f(t::Float64)::Float64
 #       his(t)/t
 #   end
    funcomp(1.0, (t)->his(t)/t, his.range)
end

fout=Operator(f0)

#axpby!(0.1,f0,10.0,fout)
#fout.f[1](1.1)
#fout.f[2](1.1)
#
##TODO: ez miért van?!?!
#fout=Operator(f0)
#axpby!(10.,fout,0.0,f0)
#f0.f[2](1.1)
#f0.f[1](1.1)
#
t=fout.range[1]:0.01:fout.range[2]
plot(t,fout.(t))
for _ in 1:5
    fout=Operator(fout)
    plot!(t,fout.(t))
end
plot!()
musFun = getindex(schursolve(Operator, f0, 3, :LM,
        KrylovKit.Arnoldi(krylovdim=10, tol=1e-1, verbosity=3)), [3, 2, 1, 4]);


#------------- eigen functions -----------------
NN = 1#1:2:7
plot(t, getindex.([musFun[2][NN](tt) for tt in t], 1), width=2, linestyle=:dash,lab=abs(musFun[1][NN]))
for NN in 2: length(musFun[2])
    plot!(t, getindex.([musFun[2][NN](tt) for tt in t], 1), width=2, linestyle=:dash,lab=abs(musFun[1][NN]))
end
plot!()
#------------- eigen functions -----------------
NN = 7#1:2:7
plot(t, getindex.([musFun[2][NN](tt)*real(musFun[1][NN]) for tt in t], 1), width=2, linestyle=:dash,lab="shape")
#Mapped:
fmapped=Operator(musFun[2][NN])
plot!(t, getindex.([fmapped(tt).*real(musFun[1][NN]) for tt in t], 1), width=2, linestyle=:dash,lab="mapped")


