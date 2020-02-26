include("jsr.jl")

A1 = [1 0; 1 0]
B1 = [0 1; 0 1]

v = Matrix{Float64}[]
push!(v, A1)
push!(v, B1)

lb, ub = jsr(v)

println("[", lb, ", ", ub, "]")


