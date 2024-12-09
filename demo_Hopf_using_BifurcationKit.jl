using BifurcationKit, Plots

falten(s) = [getindex.(s0_start, 1)..., getindex.(s0_start, 2)...]

wrap(v) = [SA[v[i], v[i+100]] for i in 1:100]


function F(x, p)

    a0 = wrap(aa)

    v0 = LinMap(dp_0_cb, s0; p=p)[1]
    return falten(v0)
end

b_H_start = 0.1
p = (ζ, δ, b_H_start, τ, μ)
p_start = p
xstar = falten(s0_start)

F(xstar, p_start)


prob = BifurcationProblem(F, xstar, p_start, 3;    record_from_solution = (x,p; k...) -> maximum(abs.(x[1:100]))))

br = continuation(prob, PALC(), ContinuationPar(p_min=b_H_start, p_max=1.5))
scene = plot(br)