
using Plots
using StaticArrays
using ArraysOfArrays
using Parameters

include("./sice_kinematics.jl")
include("./sice_robot_dynamics.jl")

using .SiceKinematics
using .SiceDynamics

function split_vec_of_arrays(u)
    vec.(u) |>
    x ->
    VectorOfSimilarVectors(x).data |>
    transpose |>
    VectorOfSimilarVectors
end


"""実験データ"""
@with_kw struct Data{T}
    t::Vector{T}
    q::Vector{Vector{T}}
    dq::Vector{Vector{T}}
    ddq::Vector{Vector{T}}
    desired_ddq::Vector{Vector{T}}
    u::Vector{Vector{T}}
    x1::Vector{Vector{T}}
    x2::Vector{Vector{T}}
    x3::Vector{Vector{T}}
    x4::Vector{Vector{T}}
    dx1::Vector{Vector{T}}
    dx2::Vector{Vector{T}}
    dx3::Vector{Vector{T}}
    dx4::Vector{Vector{T}}
    error::Vector{T}
    dit_to_obs::Vector{T}
    jl::Vector{BitVector}
    F_distur::Vector{Vector{T}}  # トルクに入る外乱
    Fc::Vector{Vector{T}}  # 対象に与えた力
end


"""全部実行"""
function run_simulation()

    TIME_SPAN = 10.0
    Δt = 0.01

    t = range(0.0, TIME_SPAN, step=Δt)

    data = Data(
        t = t,
        q = Vector{Vector{T}}(undef, length(t)),
        dq = Vector{Vector{T}}(undef, length(t)),
        ddq = Vector{Vector{T}}(undef, length(t)),
        desired_ddq = Vector{Vector{T}}(undef, length(t)),
        u = Vector{Vector{T}}(undef, length(t)),
        x1 = Vector{Vector{T}}(undef, length(t)),
        x2 = Vector{Vector{T}}(undef, length(t)),
        x3 = Vector{Vector{T}}(undef, length(t)),
        x4 = Vector{Vector{T}}(undef, length(t)),
        dx1 = Vector{Vector{T}}(undef, length(t)),
        dx2 = Vector{Vector{T}}(undef, length(t)),
        dx3 = Vector{Vector{T}}(undef, length(t)),
        dx4 = Vector{Vector{T}}(undef, length(t)),
        error = Vector{T}(undef, length(t)),
        dis_to_obs = Vector{T}(undef, length(t)),
        jl = Vector{BitVector}(undef, length(t)),
        F_distur =  Vector{Vector{T}}(undef, length(t)),
        Fc =  Vector{Vector{T}}(undef, length(t)),
    )

    # 初期値
    q₀ = [0.0, 0.0, 0.0, 0.0]
    dq₀ = zero(q)
    ddq₀ = zero(q)
    desired_ddq₀ = zero(q)
    u = zero(q)
    F_distur₀ = zero(q)

    # 目標値
    xd = [1.0, 1.0]

    # 障害物
    xo = [2.0, 2.0]

    # 初期値代入

    xi = [x1, x2, x3, x4]


    xs = [[0.0, 0.0]]
    for x in xi
        push!(xs, x(q))
    end

    x, y = split_vec_of_arrays(xs)
    plot(x, y)

    scatter!([xd[1]], [xd[2]])




end




@time run_simulation()