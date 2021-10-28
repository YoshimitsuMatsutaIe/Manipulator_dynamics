using Symbolics

@variables q[1:4]

T(q) = [
    cos(q) -sin(q)
    sin(q) cos(q)
]

r = [
    1
    2
]

# T_all(q) = T(q[1]) * T(q[2])

# x = T_all(q)
T_(q) = T(q[1])*r
x = T_(q)
J = Symbolics.jacobian(x, [q[1], q[2], q[3], q[4]])
println(J)