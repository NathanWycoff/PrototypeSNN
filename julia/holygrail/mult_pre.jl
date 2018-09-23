#  julia/forward_net.jl Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 09.10.2018

# Can we do more than 1 presynaptic action potential?

using JuMP, Gurobi

m = Model(solver=GurobiSolver())

# Params
L = 10
tau = 1
tf1 = 1
tf2 = 2
nu = 1
td = 2.5

# DV
@variable(m, b11, Bin)
@variable(m, b12, Bin)
@variable(m, b13, Bin)
@variable(m, b21, Bin)
@variable(m, b22, Bin)
@variable(m, b23, Bin)

@variable(m, ta)
@variable(m, mu)

# obj
@objective(m, Min, (ta - td)^2)

# Constr
# Enforce binary conditions
@constraint(m, -L * b11 <= ta - tf1)
@constraint(m, L * (b12 + b13) >= ta - tf1)
@constraint(m, L * b13 >= ta - tf1 - tau)
@constraint(m, -L * b21 <= ta - tf2)
@constraint(m, L * (b22 + b23) >= ta - tf2)
@constraint(m, L * b23 >= ta - tf2 - tau)

@constraint(m, b11 + b12 + b13 == 1)
@constraint(m, b21 + b22 + b23 == 1)

@constraint(m, b12 * (ta - tf1) + b13 * tau + b22 * (ta - tf2) + b23 * tau == mu * nu)

print(m)

# Solve with Gurobi
status = solve(m)

# Solution
println("Objective value: ", getobjectivevalue(m))
println("ta = ", getvalue(ta))
println("w = ", 1/getvalue(mu))
println("b11 = ", getvalue(b11))
println("b12 = ", getvalue(b12))
println("b13 = ", getvalue(b13))
println("b21 = ", getvalue(b21))
println("b22 = ", getvalue(b22))
println("b23 = ", getvalue(b23))
