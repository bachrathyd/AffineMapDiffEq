#Demo: continuouse delay and sampled delay
5 + 5
using DDE_mapping

using Plots
plotly()
using StaticArrays
using DifferentialEquations

using LinearAlgebra

using Revise
using MDBM


# Governing equation

function mixedDelay_oscill(u, h, p, t)
    # Parameters
    a, b, c, τ, dh, T = p
    #External forcing
    F = 0.0#0.1 * (cos(2pi * t / T) .^ 10)
    # Components of the delayed differential equation
    dx = u[2]
    #ddx = -(δ + ϵ * cos(2pi * t / T)) * u[1] - 2 * ζ * u[2] + b * h(p, t - τ)[1] + F
    if typeof(dh) == Float64
        tau_fix = mod(t, dh) + dh
    else
        tau_fix = mod(t, dh.value) + dh.value
    end

    ddx = -(a) * u[1] + b * h(p, t - τ)[1] + c * h(p, t - tau_fix)[1] + F
    # Update the derivative vector
    @MArray[dx, ddx]
    #SA[dx, ddx]
end
Base.:+(a::SVector, b::Bool) = a .+ b
Base.:+(a::SVector, b::Float64) = a .+ b #TODO: where to put this?


## parameters 
a = 0.6#
b = -0.2
c = -0.05
c = 0.1
τ = 2pi
dh = τ / sqrt(2)
T = dh * 5#Float64(dh) #31.4159# 2*dh#2pi * 5
p = a, b, c, τ, dh, T

#t=LinRange(0,2,5000)
#plot(t,mod.(t,Ref(dh)).+dh)
# test simulation ---------------
#initial condition
u0 = @MArray [1.0, 0.0]
#u0 = SA[1.0, 0.0]
#history function
h(p, t) = @MArray [0.0; 0.0]
#h(p, t) = SA[0.0; 0.0]
probMixtau = DDEProblem(mixedDelay_oscill, u0, h, (0.0, T * 100), p)#; constant_lags=[τ]

#Parameters for the solver as a Dict (it is necessary to collect it for later use)
#Solver_args = Dict(:alg => MethodOfSteps(BS3()), :callback => cb, :adaptive => true, :dt => 0.01, :verbose => false, :reltol => 1e-9)#, save_everystep=false)#abstol,reltol)
Solver_args = Dict(:alg => MethodOfSteps(Tsit5()), :verbose => false, :reltol => 1e-7)#
#Solver_args = Dict()#
sol = solve(probMixtau; Solver_args...)#abstol,reltol


plot(sol)

#last period of the long simulation:
t_select_period = 0.0:0.01:T
t_select_delay = eriod = 0.0:0.001:τ
sol_period = sol(sol.t[end] .- t_select_period)
sol_delay = sol(sol.t[end] .- t_select_delay)

#plot the state phase of the last segment
plot(sol_period.t, getindex.(sol_period.u, 1))
plot!(sol_period.t, getindex.(sol_period.u, 2))
plot!(sol_delay.t, getindex.(sol_delay.u, 1))
plot!(sol_delay.t, getindex.(sol_delay.u, 2))
#plot the phase space (u - du)
plot(getindex.(sol_delay.u, 1), getindex.(sol_delay.u, 2))
plot!(getindex.(sol_period.u, 1), getindex.(sol_period.u, 2))


# ---------------- Affine mapping ---------------------------
using KrylovKit
Neig = 8#number of required eigen values
Krylov_arg = (Neig, :LM, KrylovKit.Arnoldi(tol=1e-32, krylovdim=18 + 5, verbosity=0));

τmax = maximum([τ, 2 * dh]) * 1.1 #maximal timedelay in the mapping
Nstep = 5 # discretization number of the mapping
Timeperiod = 5 * dh# dh*2#copy(dh) # timeperiod of the mapping

#Creating the problem
probMixtau = DDEProblem(mixedDelay_oscill, u0, h, (0.0, Timeperiod), p)#; constant_lags=[τ]

dp_mix_delay = dynamic_problemSampled(probMixtau, Solver_args, τmax,
    Timeperiod; Historyresolution=Nstep,
    zerofixpont=true, affineinteration=8,
    Krylov_arg=Krylov_arg)


#solvig the problem (mu: Floquet multiplier, saff: discretized foxed point (in [-τmax,0]), sol: the corresponding periodic solution of the fixpoint)
@time mu, saff, sol0 = affine(dp_mix_delay; p=(a, b, c, τ, dh, T));



dp = dp_mix_delay

hint(p, t) = @MArray [0.1; 1.0]
sol = solve(remake(dp.Problem; u0=hint(p, 0.0), h=hint, p=(a, b, c, τ, dh, T)); dp.alg...)

# Comparing the solutions:
plot(sol_period[1, :], sol_period[2, :], marker=:circle, markersize=6, lab="")
plot!(getindex.(saff, 1), getindex.(saff, 2), marker=:circle, markersize=4, lab="")#
plot!(sol0[1, :], sol0[2, :], marker=:circle, lw=0, lab="")#marker=:cross,markersize=2)#
#


# Plotting the Floquet multipliers
#in log scale
plot(log.(abs.(mu[1])))
#in complex plane
scatter(mu[1])
plot!(sin.(0:0.01:2pi), cos.(0:0.01:2pi))








## ----------------------- creating stability chart -------------------------

Krylov_arg = (3, :LM, KrylovKit.Arnoldi(tol=1e-12, krylovdim=12, verbosity=0));

Nstep = 20
dp_mix_delay = dynamic_problemSampled(probMixtau, Solver_args, τmax,
    Timeperiod; Historyresolution=Nstep,
    zerofixpont=true, affineinteration=0,
    Krylov_arg=Krylov_arg)

# ------------------ Peak-2-Peak ampNorm color, and spectral radius - BRUTE FORCE in the plane of δ-ϵ  ----------------

av = -0.2:0.01:1.0 # initial grid in x direction
bv = -0.5:0.01:0.5 # initial grid in y direction
Aaff = zeros(size(av, 1), size(bv, 1))
Spek_aff = zeros(size(av, 1), size(bv, 1))

@time Threads.@threads for j in 1:size(bv, 1)
    @inbounds b = bv[j]
    # @show j
    Threads.@threads for i in 1:size(av, 1)
        @inbounds a = av[i]
        muaff, s0aff, sol0 = affine(dp_mix_delay; p=(a, b - c, c, τ, dh, T))
        Aaff[i, j] = norm(getindex.(s0aff, 1)) # norm of the motion
        # Aaff[i,j]= maximum(abs.(getindex.(s0aff,1))) #maximum amplitude of the orbit
        Spek_aff[i, j] = maximum(abs.(muaff[1]))
    end
end


#Plotting the maximal amplitud on the stable domain only
Spek_affsat = deepcopy(Spek_aff);
Spek_affsat[Spek_affsat.>1.0] .= 0.0;
#heatmap(bv, av, log.(Spek_affsat))
heatmap(av, bv, log.(Spek_affsat'))


#------------------ Stability boundary - MDBM -----------------

ax1 = Axis(-1.0:0.1:1.0, "a") # initial grid in x direction
ax2 = Axis(-0.5:0.25:0.5, "b") # initial grid in y direction
function fooDelay(a, b)
    ABSmuMax = spectralradius(dp_mix_delay; p=(a, b - c, c, τ, dh, T))
    #return ABSmuMax - 1.0
    return log(ABSmuMax)
end
mymdbm = MDBM_Problem(fooDelay, [ax1, ax2], memoization=true)
@time MDBM.solve!(mymdbm, 4)

#points where the function foo was evaluated
x_sol, y_sol = getinterpolatedsolution(mymdbm)
#scatter(x_eval,y_eval,markersize=1)
#scatter!(x_sol,y_sol,markersize=2)
scatter!(x_sol, y_sol, markersize=1)

# #------------------ Stability boundary - MDBM -----------------
# 
# ax1 = Axis(-1.0:0.5:1.0, "a") # initial grid in x direction
# ax2 = Axis(-0.5:0.2:0.5, "b") # initial grid in y direction
# ax3 = Axis(-0.3:1:pi+0.5, "ω") # initial grid in y direction
# function fooDelay(a, b, ω)
# 
#     c = -0.05
#     τ = 2pi
#     dh = τ
#     # dh = τ / sqrt(2) -don't work
#     μ = exp(1im * ω)
#     #(2*a*μ^2-(c+2*b)*μ-c)  *cosh(h*sqrt(b/μ-a)) = a*μ ^3-b*μ^2+(a-c)*μ-c-b
#     D = (2 * a * μ^2 - (c + 2 * b) * μ - c) * cosh(dh * sqrt(b / μ - a)) - (a * μ^3 - b * μ^2 + (a - c) * μ - c - b)
#     return imag(D), real(D)
# end
# mymdbm = MDBM_Problem(fooDelay, [ax1, ax2, ax3], memoization=true)
# @time MDBM.solve!(mymdbm, 5)
# #points where the function foo was evaluated
# xyz_sol = getinterpolatedsolution(mymdbm)
# 
# scatter!(xyz_sol[1], xyz_sol[2], markersize=2, size=(1200, 800))
# 

##--------------------------

##--------------------------

#------------------ Stability boundary - MDBM -----------------

ax1 = Axis(-1.0:0.5:2.7, "a") # initial grid in x direction
ax2 = Axis(-1.5:0.5:1.5, "b") # initial grid in y direction
ax3 = Axis(-0.0:0.5:2.0, "c") # initial grid in y direction
#ax3 = Axis(-2.0:0.5:0.0, "c") # initial grid in y direction

ax1 = Axis(-0.5:0.5:1.0, "a") # initial grid in x direction
ax2 = Axis(-0.5:0.25:0.5, "b") # initial grid in y direction
ax3 = Axis(-0.6:0.2:0.6, "c") # initial grid in y direction
#ax3 = Axis(-0.6:0.1:-0.0, "c") # initial grid in y direction
#function fooDelay(a, b, c)
function fooDelay(a::Float64, b::Float64, c::Float64)::Float64
    ABSmuMax = spectralradius(dp_mix_delay; p=(a, b, c, τ, dh, T))
    #return ABSmuMax - 1.0
    return log(ABSmuMax)
end
mymdbm3D = MDBM_Problem(fooDelay, [ax1, ax2, ax3], memoization=true)
@time MDBM.solve!(mymdbm3D, 1, verbosity=1)

#points where the function foo was evaluated
xyz_sol = getinterpolatedsolution(mymdbm3D)
#scatter(x_eval,y_eval,markersize=1)
#scatter!(x_sol,y_sol,markersize=2)
scatter(xyz_sol..., markersize=1)

#scatter!(x_sol, y_sol ,c .+ 0.0 .* y_sol, markersize=1)



#--------------------------- Sub-cube interpolation----------------
#calcuatin the sub-cubes interpolations stored in the mymdbm.ncubes[i].posinterp
MDBM.interpsubcubesolution!(mymdbm3D)
#extracting the resutls to from the 
path2points = extract_paths(mymdbm3D)

#extracting the unique points and plotting
puniq = unique(collect(Iterators.flatten(Iterators.flatten(path2points))))
scatter!( getindex.(puniq, 1), getindex.(puniq, 2), getindex.(puniq, 3), markersize=3, color=:green,label = "subface - solution")



#exctracing the simplexes for each ncube
flatened_path2points = collect(Iterators.flatten(path2points))
#eliminating the points with less than 2 points (caused by fininte precision)
true_truflatened_path2points = flatened_path2points[length.(flatened_path2points).==3]


#plotting the lines between the points
using GLMakie
GLMakie.closeall()
GLMakie.activate!(; title="3 parameters, codimension 1")
f = Figure(size=(1300, 800))
lights = [
    DirectionalLight(RGBf(0, 0, 0.7), Vec3f(-1, -1, 0)),
    DirectionalLight(RGBf(0.7, 0.2, 0), Vec3f(-1, 1, -1)),
    DirectionalLight(RGBf(0.7, 0.7, 0.7), Vec3f(1, -1, -1))
]
ax2 = LScene(f[1, 1], scenekw=(lights=lights,))

n_faces = reshape(1:(3*length(true_truflatened_path2points)), (3, length(true_truflatened_path2points)))'
vertices_mat = hcat(Iterators.flatten(true_truflatened_path2points)...)
GLMakie.mesh!(ax2, vertices_mat, n_faces, color=:steelblue, alpha=0.5,
    label="subface - local simplex")

axislegend(ax2)

display(GLMakie.Screen(), f)


##


using LinearAlgebra
function export_mesh_to_stl(filename::String, vertices_mat::AbstractMatrix, n_faces::AbstractMatrix)
    open(filename, "w") do io
        println(io, "solid mesh")

        for j in 1:size(n_faces, 2)
            # Get indices of triangle vertices (assuming 0-based indices from GLMakie)
            i1, i2, i3 = n_faces[:, j]

            # Get triangle vertices
            v1 = vertices_mat[:, i1]
            v2 = vertices_mat[:, i2]
            v3 = vertices_mat[:, i3]

            # Compute normal (normalized cross product)
            n = cross(v2 - v1, v3 - v1)
            norm_n = norm(n)
            n = norm_n == 0 ? [0.0, 0.0, 0.0] : n / norm_n

            # Write triangle
            println(io, "  facet normal $(n[1]) $(n[2]) $(n[3])")
            println(io, "    outer loop")
            println(io, "      vertex $(v1[1]) $(v1[2]) $(v1[3])")
            println(io, "      vertex $(v2[1]) $(v2[2]) $(v2[3])")
            println(io, "      vertex $(v3[1]) $(v3[2]) $(v3[3])")
            println(io, "    endloop")
            println(io, "  endfacet")
        end

        println(io, "endsolid mesh")
    end
end


using LinearAlgebra

function export_mesh_to_binary_stl(filename::String, vertices::AbstractMatrix, faces::AbstractMatrix)
    n_tri = size(faces, 2)
    open(filename, "w") do io
        # 80-byte header
        hdr = zeros(UInt8, 80)
        txt = "Binary STL from Julia"
        hdr[1:length(txt)] = codeunits(txt)
        write(io, hdr)

        # write triangle count (UInt32 little-endian)
        write(io, UInt32(n_tri))

        for j in 1:n_tri
            i1, i2, i3 = faces[:, j]
            v1, v2, v3 = vertices[:, i1], vertices[:, i2], vertices[:, i3]
            n = cross(v2 - v1, v3 - v1)
            if norm(n) > 0
                n ./= norm(n)
            else
                n .= 0.0
            end

            # write normal + verts as Float32 LE
            write(io,
                Float32(n[1]), Float32(n[2]), Float32(n[3]),
                Float32(v1[1]), Float32(v1[2]), Float32(v1[3]),
                Float32(v2[1]), Float32(v2[2]), Float32(v2[3]),
                Float32(v3[1]), Float32(v3[2]), Float32(v3[3]),
                UInt16(0)           # attribute byte count
            )
        end
    end
end

export_mesh_to_stl("output_h_sqrt2_tau.stl", vertices_mat, n_faces')
export_mesh_to_binary_stl("output_binary_h_sqrt2_tau.stl", vertices_mat, n_faces')