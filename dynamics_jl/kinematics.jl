using CPUTime
using Plots

include("./lagrange.jl")


joint_positions = Matrix(undef, 3, n)

for i in 1:n
    t = T(1, i)
    joint_positions[:, i] = t[1:3, 4]
end


plot(
    joint_positions[1, :], joint_positions[2, :], joint_positions[3, :],
    aspect_ratio = 1,
    marker=:circle, label = "joints",
)
