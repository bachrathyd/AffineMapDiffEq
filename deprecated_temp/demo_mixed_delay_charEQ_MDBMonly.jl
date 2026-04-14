#Demo: continuouse delay and sampled delay
5 + 5
using Plots
plotly()
using LinearAlgebra
using MDBM

ax1 = Axis(-1.0:0.5:1.0, "a") # initial grid in x direction
ax2 = Axis(-0.5:0.2:0.5, "b") # initial grid in y direction
ax3 = Axis(-0.3:0.2:pi+0.3, "ω") # initial grid in y direction
function fooDelay(a::Float64, b::Float64, ω::Float64)::Tuple{Float64, Float64}
    c = -0.05
    τ = 2pi
    dh = τ
    μ = exp(1im * ω)
    #(2*a*μ^2-(c+2*b)*μ-c)  *cosh(h*sqrt(b/μ-a)) = a*μ ^3-b*μ^2+(a-c)*μ-c-b
    D = (2*a*μ^2-(c+2*b)*μ-c)  *cosh(dh*sqrt(b/μ-a)) -( a*μ ^3-b*μ^2+(a-c)*μ-c-b)
    (imag(D), real(D))
end

mymdbm = MDBM_Problem(fooDelay, [ax1,ax2,ax3],memoization=false)
doThreadprecomp=false
@time solve!(mymdbm, 1)
@time solve!(mymdbm, 1)
@time solve!(mymdbm, 1)
@time solve!(mymdbm, 1)

#points where the function foo was evaluated
xyz_sol = getinterpolatedsolution(mymdbm)

aaa = scatter(xyz_sol..., markersize=2)
display(aaa)

bbb=scatter(xyz_sol[1:2]..., markersize=2)

# ## ----------------------------------------  
# 
# #Demo: continuouse delay and sampled delay
# 5 + 5
# using Plots
# plotly()
# using LinearAlgebra
# using MDBM
# 
# ax1 = Axis(-1.0:0.5:1.0, "a") # initial grid in x direction
# ax2 = Axis(-0.5:0.2:0.5, "b") # initial grid in y direction
# ax3 = Axis(-0.1:1:2pi+0.8, "ω") # initial grid in y direction
# function fooDelay(a::Float64, b::Float64, ω::Float64)::Tuple{Float64,Float64}
#     c = -0.05
#     τ = 2pi
#     dh = τ
#     μ = exp(1im * ω)
#     #(2*a*μ^2-(c+2*b)*μ-c)  *cosh(h*sqrt(b/μ-a)) = a*μ ^3-b*μ^2+(a-c)*μ-c-b
#     D = (2 * a * μ^2 - (c + 2 * b) * μ - c) * cosh(dh * sqrt(b / μ - a)) - (a * μ^3 - b * μ^2 + (a - c) * μ - c - b)
#     return imag(D), real(D)
# end
# #@code_warntype fooDelay(1.1,2.2,3.3)
# mymdbm = MDBM_Problem(fooDelay, [ax1, ax2, ax3], memoization=true)
# @time MDBM.solve!(mymdbm, 1)
# #points where the function foo was evaluated
# xyz_sol = getinterpolatedsolution(mymdbm)
# 
# #scatter+(xyz_sol[1], xyz_sol[2], markersize=2)
# #5+5
# 
# 
# aaa = scatter(xyz_sol..., markersize=2)
# display(aaa)