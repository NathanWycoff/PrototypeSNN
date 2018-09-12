#  julia/neur_playground.jl Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 09.09.2018

# Can we do it?

using JuMP, Gurobi

m = Model(solver=GurobiSolver())

# Params
L = 10
tau = 1
tf = 1
nu = 1
td = 1.1

# DV
@variable(m, b1, Bin)
@variable(m, b2, Bin)
@variable(m, b3, Bin)

@variable(m, ta)
@variable(m, mu)

# obj
@objective(m, Min, (ta - td)^2)

# Constr
@constraint(m, -L * b1 <= ta - tf)
@constraint(m, L * b2 >= ta - tf)
@constraint(m, L * b3 >= ta - tf - tau)
@constraint(m, b1 + b2 + b3 == 1)

@constraint(m, b2 * (ta - tf) + b3 * tau == mu * nu)

print(m)

# Solve with Gurobi
status = solve(m)

# Solution
println("Objective value: ", getobjectivevalue(m))
println("ta = ", getvalue(ta))
println("w = ", 1/getvalue(mu))
println("b1 = ", getvalue(b1))
println("b2 = ", getvalue(b2))
println("b3 = ", getvalue(b3))
