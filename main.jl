include("jsr.jl")

#A1 = [1 0; 1 0]
#B1 = [0 1; 0 1]

A1 = [0.8147 0.1270; 0.9058 0.9134]
B1 = [0.6324 0.2785; 0.0975 0.5469]
# JSR: bounds on the jsr: [1.20676933149348, 1.20676945215481]

A2 = [-0.2565 -0.1500; 0.4293 -0.3034]
B2 = [-0.2489 -0.0267; 0.1160 -0.1483]
# JSR: bounds on the jsr: [0.377110401321558, 0.37711043903225]

v = Matrix{Float64}[]
push!(v, A1)
push!(v, B1)
lb, ub = jsr(v, verbose = true)

v = Matrix{Float64}[]
push!(v, A2)
push!(v, B2)
lb, ub = jsr(v, verbose = true)

