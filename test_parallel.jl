

# test - Parallel 
t=0:0.01:2*pi
#x0=t
x0=(mod.(t,2*pi/2)) .^ 1.0
x0[x0 .> 1.0] .= 1.0
x0[x0 .< 1.0] .= 0.0
plot(t,x0)

for kk in 0:200
   ## kk+=1
    ##kk=1000
xout=remove_parallel_components(x0, [sin.(k .* t) for k in 0:kk])
xout=remove_parallel_components(xout, [cos.(k .* t) for k in 0:kk])
plot!(t,xout)
println(norm(xout))
end
plot!(t,xout)
