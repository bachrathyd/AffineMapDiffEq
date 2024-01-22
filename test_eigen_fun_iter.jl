

# testing fix point and eigen function interation

v0i=fix_v
#v0i=s0
ks=0
plot(sol)
for _  in 1:1
v0i=OneMap(v0i)
#v0i=remove_parallel_components(v0i, mus[2])
ks +=1
end
plot((ks*Nstep .+ (1:Nstep) ) ./ Nstep,real.(vfix_simulation))
plot!((ks*Nstep .+ (1:Nstep)) ./ Nstep ,real.(v0i))
#plot!(imag.(v0i))

II=1
plot(real.(mus[2][II]))
plot!(imag.(mus[2][II]))

A1=mus[2][II]
A1=OneMap(A1)
plot!(real.(A1 ./ mus[1][II]))
plot!(imag.(A1 ./ mus[1][II]))

plot(sol.t,sol.u)
for tt= -3:0.2:5
    xx= interpolate_complex_on_grid(x0,  t[1], t[2]-t[1],tt)
scatter!([tt],[xx])
end
scatter!(0,0)

