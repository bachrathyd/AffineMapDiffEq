#Longitudinal_vib_PDF with feedback from the boundaries

# WORKING verstion!!!!!!!!!!!!!!!!!

5 + 5
using ModelingToolkit, MethodOfLines, OrdinaryDiffEq, DifferentialEquations, DomainSets

using Plots



# Define constants (replace with your specific values)
c=0.1#sound speed


η0 = 0.02#5
η(x) = η0#damping

L = 10.0 # Length of the beam
Tend = 200.5 # end time

@parameters x t
@variables w(..), F(..)
Dt = Differential(t)
Dtt = Differential(t)^2
Dx = Differential(x)
Dxx = Differential(x)^2
Dxxx = Differential(x)^3
K=-.8900
tau=0.01
ω = 10.0*2*pi
q(x, t) = 0.0 * cos(ω * t)

w0(x) = 0.01 * exp(-((x - L / 5) / L * 100/4)^2)
w0(x, t) = 1.0*w0(x)
#BC_time(t) = 0.01 * exp(-((t -0.02) / 0.005)^2)
#BC_time(x, t)  = BC_time(t)
#display(plot(w0.(LinRange(0, L, 500))))
#display(plot(BC_time.(LinRange(0, Tend, 500))))

# Define spatial domain and grid
x_min = t_min = 0.0
x_max = L
t_max = Tend

#eq = [ρA * Dtt(w(x, t)) + η(x)*Dt(w(x, t)) + Dxx(EI*(Dxx(w(x, t))))~ + q(x, t) ]
eq = [ Dtt(w(x, t))+ η0*Dt(w(x, t))  - c * Dxx(w(x, t)) ~ 0, F(t) ~ Dx(w(0, t))]



domains = [x ∈ Interval(x_min, x_max),
    t ∈ Interval(t_min, t_max)]


# #Boundary conditions #Free-Free
# bcs = [ w(x, 0) ~ w0(x, t), #Initail Condition: position
# Dt(w(x, 0)) ~ 0.0, #Initail Condition: velocity
# Dx(w(0, t)) ~  0.0, # Boundary Condition: left
# Dx(w(x_max, t)) ~  0.0,]  # Boundary Condition: right

# #Boundary conditions #Clamped-Free 
# bcs = [ w(x, 0) ~ w0(x, t), #Initail Condition: position
#        Dt(w(x, 0)) ~ 0.0, #Initail Condition: velocity
#        w(0, t) ~  0.0, # Boundary Condition: left
#        Dx(w(x_max, t)) ~  0.0,]  # Boundary Condition: right

# # Boundary conditions - #Clamped-Clamped
# bcs = [w(x, 0) ~ w0(x, t), #Initail Condition: position
#    Dt(w(x, 0)) ~ 0.0,#Initail Condition: velocity
#    w(0, t) ~ 0.0,# Boundary Condition: left
#    w(x_max, t) ~ 0.0,]# Boundary Condition: right
#    #Dx(w(x_max, t)) ~ Dx(w(x_max, t)),]# Boundary Condition: right
#    #Dx(w(x_max, t)) ~ Dx(w(x_max, t-tau)),]# Boundary Condition: right


# #Boundary conditions #Free - Clamped
# bcs = [ w(x, 0) ~ w0(x, t), #Initail Condition: position
# Dt(w(x, 0)) ~ 0.0, #Initail Condition: velocity
# Dx(w(0, t)) ~  0.0, # Boundary Condition: left
# w(x_max, t) ~  0.0,]  # Boundary Condition: right



 #Boundary conditions #Delayed cotrolled - Clamped
 bcs = [ w(x, 0) ~ w0(x, t), #Initail Condition: position
 Dt(w(x, 0)) ~ 0.0, #Initail Condition: velocity
 w(0, t) ~  0.0,# Boundary Condition: left
 Dx(w(x_max, t)) ~ K * F(t) , 
 F(0) ~ Dx(w(x_max, 0))
]  # Boundary Condition: right



@named pdesys = PDESystem(eq, bcs, domains, [x, t], [w(x, t),F(t)])


N = 100  # Number of spatial grid points
order = 4  # Order of the spatial discretization
discretization = MOLFiniteDifference([x => N], t, approx_order=order)


println("Discretization:")
@time prob = discretize(pdesys, discretization)



println("Solve:")
t_save_plot=LinRange(0, Tend, 3500)
@time sol = solve(prob, Tsit5(), reltol=1e-3, dtmax=1e-2,saveat=t_save_plot)

discrete_x = sol[x]
discrete_t = sol[t]

solw = sol[w(x, t)]
solF = sol[F(t)]
@show size(solw)


#plotly()
gr()
N_plot_filter=1
solw[solw .> 0.03] .= 0.0
display(heatmap(discrete_t[1:N_plot_filter:end],discrete_x,solw[:,1:N_plot_filter:end]))
plot(discrete_t,solF)
##

# Make an animation0
anim = @animate for k in 1:10:(length(discrete_t)÷1)
    plot(discrete_x,solw[1:end, k], title="$(discrete_t[k])",ylim=(-0.01,0.01)) # 2:end since end = 1, periodic condition
end
display(gif(anim, "Beam_vibration_1.gif", fps = 15))

println("-------------------------------- DONE --------------------------------")
### 
### ##
### @time sol_no_wrap = solve(prob, Tsit5(), reltol=1e-3, dtmax=1e-2; wrap = Val(false))
### plot(sol_no_wrap(0.1)[501:end])
### anim = @animate for tloc in LinRange(0, Tend, 500)
###     plot(sol_no_wrap(tloc)[500:end], title="$tloc",ylim=(-0.01,0.01)) # 2:end since end = 1, periodic condition
### end
### display(gif(anim, "Beam_vibration_time.gif", fps = 25))
### ##
### 

