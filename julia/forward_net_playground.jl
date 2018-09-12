#  /Users/nathw95/Documents/Research/spiking_nets/prototype/julia/forward_net_playground.jl Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 09.12.2018

# Can we do more than 1 presynaptic action potential?

using JuMP, Gurobi

########################### 
########################### 
########################### Try with 2 presynaptic fires
m = Model(solver=GurobiSolver())

# Params
L = 10
tau = 1
tf1 = 1
tf2 = 2
nu = 1
td = 2.5
w = 0.7

# DV
# bij belongs to the ith firing time
@variable(m, b11, Bin)
@variable(m, b12, Bin)
@variable(m, b13, Bin)
@variable(m, b21, Bin)
@variable(m, b22, Bin)
@variable(m, b23, Bin)

@variable(m, ta >= 0)

# obj
#@objective(m, Min, (ta - td)^2)
@objective(m, Min, 0)

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

@constraint(m, w*(b12 * (ta - tf1) + b13 * tau + b22 * (ta - tf2) + b23 * tau ) == nu)

print(m)

# Solve with Gurobi
status = solve(m)

# Solution
println("Objective value: ", getobjectivevalue(m))
println("ta = ", getvalue(ta))
println("b11 = ", getvalue(b11))
println("b12 = ", getvalue(b12))
println("b13 = ", getvalue(b13))
println("b21 = ", getvalue(b21))
println("b22 = ", getvalue(b22))
println("b23 = ", getvalue(b23))

########################### 
########################### 
########################### 
########################### Try a 2 neuron system, each neuron fires once
m = Model(solver=GurobiSolver())

# Params
L = 10
tau = 1
tf1 = 1
nu = 1
td = 2.5
w1 = 1.5
w2 = 1.5

# DV
# bij belongs to the ith neuron
@variable(m, b11, Bin)
@variable(m, b12, Bin)
@variable(m, b13, Bin)
@variable(m, b21, Bin)
@variable(m, b22, Bin)
@variable(m, b23, Bin)

@variable(m, ta1 >= 0)
@variable(m, ta2 >= 0)

# obj
#@objective(m, Min, (ta - td)^2)
@objective(m, Min, 0)

# Constr
# Enforce binary conditions
@constraint(m, -L * b11 <= ta1 - tf1)
@constraint(m, L * (b12 + b13) >= ta1 - tf1)
@constraint(m, L * b13 >= ta1 - tf1 - tau)
@constraint(m, -L * b21 <= ta2 - ta1)
@constraint(m, L * (b22 + b23) >= ta2 - ta1)
@constraint(m, L * b23 >= ta2 - ta1 - tau)

@constraint(m, b11 + b12 + b13 == 1)
@constraint(m, b21 + b22 + b23 == 1)

@constraint(m, w1*(b12 * (ta1 - tf1) + b13 * tau) == nu)
@constraint(m, w2*(b22 * (ta2 - ta1) + b23 * tau) == nu)

print(m)

# Solve with Gurobi
status = solve(m)

# Solution
println("Objective value: ", getobjectivevalue(m))
println("ta1 = ", getvalue(ta1))
println("ta2 = ", getvalue(ta2))
