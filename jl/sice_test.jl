
using Plots
using StaticArrays
using ArraysOfArrays
using Parameters
using LinearAlgebra

include("./sice_kinematics.jl")
include("./sice_robot_dynamics.jl")
include("./rmp.jl")

using .SiceKinematics
using .SiceDynamics


const q_min = [-5/4*pi, -5/4*pi, -5/4*pi, -5/4*pi]
const q_max = [5/4*pi, 5/4*pi, 5/4*pi, 5/4*pi]



function split_vec_of_arrays(u)
    vec.(u) |>
    x ->
    VectorOfSimilarVectors(x).data |>
    transpose |>
    VectorOfSimilarVectors
end

function calc_min_dis_to_obs(xs::Vector{Vector{T}}, xo::Vector{T}) where T

    d = Vector{T}(undef, length(xs))
    for (i, x) in enumerate(xs)
        d[i] = norm(x - xo)
    end
    
    return minimum(d)
end

"""ジョイント制限を守れてるかチェック"""
function check_JointLimitation(q::Vector{T}) where T
    return (q .< q_min) .| (q_max .< q)
end

"""実験データ"""
struct Data{T}
    t::StepRangeLen{T}
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
    J1::Vector{Matrix{T}}
    J2::Vector{Matrix{T}}
    J3::Vector{Matrix{T}}
    J4::Vector{Matrix{T}}
    dJ1::Vector{Matrix{T}}
    dJ2::Vector{Matrix{T}}
    dJ3::Vector{Matrix{T}}
    dJ4::Vector{Matrix{T}}
    f0::Vector{Vector{T}}
    f1::Vector{Vector{T}}
    f2::Vector{Vector{T}}
    f3::Vector{Vector{T}}
    f4::Vector{Vector{T}}
    M0::Vector{Matrix{T}}
    M1::Vector{Matrix{T}}
    M2::Vector{Matrix{T}}
    M3::Vector{Matrix{T}}
    M4::Vector{Matrix{T}}
    error::Vector{T}
    min_dit_to_obs::Vector{T}
    jl::Vector{BitVector}
    F_distur::Vector{Vector{T}}  # トルクに入る外乱
    Fc::Vector{Vector{T}}  # 対象に与えた力
end



"""全部実行"""
function run_simulation(TIME_SPAN::T=0.1, Δt::T=0.01) where T


    # rmpのパラメータ
    attractor = RMPfromGDSAttractor(
        max_speed = 1.0,
        gain = 5.0,
        f_alpha = 1.0,
        sigma_alpha = 1.0,
        sigma_gamma = 1.0,
        wu = 5.0,
        wl = 1.0,
        alpha = 1.0,
        epsilon = 0.5,
    )

    obs_avoidance = RMPfromGDSCollisionAvoidance(
        rw = 0.1,
        sigma = 1.0,
        alpha = 1.0,
    )

    jl_avoidance = RMPfromGDSJointLimitAvoidance(
        gamma_p = 1.0,
        gamma_d = 1.0,
        lambda = 1.0,
        sigma = 1.0,
    )

    # 初期値
    q₀ = zeros(T, 4)
    dq₀ = zero(q₀)
    ddq₀ = zero(q₀)
    desired_ddq₀ = zero(q₀)
    u₀ = zero(q₀)
    F_distur₀ = zero(q₀)

    # 目標値
    xd = [1.0, 1.0]

    # 障害物
    xo = [[2.0, 2.0]]
    t = range(0.0, TIME_SPAN, step=Δt)

    data = Data(
        t,
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Matrix{T}}(undef, length(t)),
        Vector{Matrix{T}}(undef, length(t)),
        Vector{Matrix{T}}(undef, length(t)),
        Vector{Matrix{T}}(undef, length(t)),
        Vector{Matrix{T}}(undef, length(t)),
        Vector{Matrix{T}}(undef, length(t)),
        Vector{Matrix{T}}(undef, length(t)),
        Vector{Matrix{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Matrix{T}}(undef, length(t)),
        Vector{Matrix{T}}(undef, length(t)),
        Vector{Matrix{T}}(undef, length(t)),
        Vector{Matrix{T}}(undef, length(t)),
        Vector{Matrix{T}}(undef, length(t)),
        Vector{T}(undef, length(t)),
        Vector{T}(undef, length(t)),
        Vector{BitVector}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
    )



    # 初期値代入
    data.q[1] = q₀
    data.dq[1] = dq₀
    data.ddq[1] = ddq₀
    data.desired_ddq[1] = desired_ddq₀
    data.u[1] = u₀

    data.x1[1], data.x2[1], data.x3[1], data.x4[1] = calc_x(q₀)
    data.J1[1], data.J2[1], data.J3[1], data.J4[1] = calc_J(q₀)
    data.dJ1[1], data.dJ2[1], data.dJ3[1], data.dJ4[1] = calc_dJ(q₀, dq₀)

    data.dx1[1] = data.J1[1] * data.dq[1]
    data.dx2[1] = data.J2[1] * data.dq[1]
    data.dx3[1] = data.J3[1] * data.dq[1]
    data.dx4[1] = data.J4[1] * data.dq[1]

    data.f0[1] = zero(q₀)
    data.f1[1] = zeros(T, 2)
    data.f2[1] = zeros(T, 2)
    data.f3[1] = zeros(T, 2)
    data.f4[1] = zeros(T, 2)

    data.M0[1] = zeros(T, 4, 4)
    data.M1[1] = zeros(T, 2, 2)
    data.M2[1] = zeros(T, 2, 2)
    data.M3[1] = zeros(T, 2, 2)
    data.M4[1] = zeros(T, 2, 2)

    data.error[1] = norm(data.x4[1] - xd)
    data.min_dit_to_obs[1] = calc_min_dis_to_obs(
        [data.x1[1], data.x2[1], data.x3[1], data.x4[1]],
        xo
    )
    data.jl[1] = check_JointLimitation(q₀)
    data.F_distur[1] = zero(u₀)
    data.Fc[1] = zero(u₀)


    # ぐるぐる回す
    for i in 1:length(data.t)-1
        # 状態を更新
        data.q[i+1] = data.q[i] .+ data.dq[i] .* Δt
        data.dq[i+1] = data.dq[i] .+ data.ddq[i] .*Δt

        data.x1[i+1], data.x2[i+1], data.x3[i+1], data.x4[i+1] = calc_x(data.q[i+1])
        data.J1[i+1], data.J2[i+1], data.J3[i+1], data.J4[i+1] = calc_J(data.q[i+1])
        data.dJ1[i+1], data.dJ2[i+1], data.dJ3[i+1], data.dJ4[i+1] = calc_dJ(data.q[i+1], data.dq[i+1])

        data.dx1[i+1] = data.J1[i+1] * data.dq[i+1]
        data.dx2[i+1] = data.J2[i+1] * data.dq[i+1]
        data.dx3[i+1] = data.J3[i+1] * data.dq[i+1]
        data.dx4[i+1] = data.J4[i+1] * data.dq[i+1]



        # rmpを計算

        # rmpをセット
        # root
        data.f0[i+1], data.M0[i+1] = get_natural(
            jl_avoidance, data.q[i+1], data.dq[i+1], q_max, q_min, q₀
        )
        data.f1[i+1] = zeros(T, 2)
        data.f2[i+1] = zeros(T, 2)
        data.f3[i+1] = zeros(T, 2)
        data.f4[i+1] = zeros(T, 2)
        data.M1[i+1] = zeros(T, 2, 2)
        data.M2[i+1] = zeros(T, 2, 2)
        data.M3[i+1] = zeros(T, 2, 2)
        data.M4[i+1] = zeros(T, 2, 2)


        for o in xo
            f, M = get_natural(obs_avoidance, data.x1[i+1], data.dx1[i+1], o)
            @. data.f1[i+1] += f
            @. data.M1[i+1] += M

            f, M = get_natural(obs_avoidance, data.x2[i+1], data.dx2[i+1], o)
            @. data.f2[i+1] += f
            @. data.M2[i+1] += M

            f, M = get_natural(obs_avoidance, data.x3[i+1], data.dx3[i+1], o)
            @. data.f3[i+1] += f
            @. data.M3[i+1] += M

            f, M = get_natural(obs_avoidance, data.x4[i+1], data.dx4[i+1], o)
            @. data.f4[i+1] += f
            @. data.M4[i+1] += M
        end

        f, M = get_natural(attractor, data.x4[i+1], data.dx4[i+1], xd)
        @. data.f4[i+1] += f
        @. data.M4[i+1] += M


        # pullback演算
        root_f = data.f0[i+1]
        root_M = data.M0[i+1]

        _rf, _rM = pullbacked_rmp(data.f1[i+1], data.M1[i+1], data.J1[i+1], data.dJ1[i+1], data.dx1[i+1])
        @. root_f += _rf
        @. root_M += _rM

        _rf, _rM = pullbacked_rmp(data.f2[i+1], data.M2[i+1], data.J2[i+1], data.dJ2[i+1], data.dx2[i+1])
        @. root_f += _rf
        @. root_M += _rM

        _rf, _rM = pullbacked_rmp(data.f3[i+1], data.M3[i+1], data.J3[i+1], data.dJ3[i+1], data.dx3[i+1])
        @. root_f += _rf
        @. root_M += _rM

        _rf, _rM = pullbacked_rmp(data.f4[i+1], data.M4[i+1], data.J4[i+1], data.dJ4[i+1], data.dx4[i+1])
        @. root_f += _rf
        @. root_M += _rM

        # resolve演算
        data.desired_ddq[i+1] = pinv(root_M) * root_f

        data.error[1] = norm(data.x4[1] - xd)
        data.min_dit_to_obs[1] = calc_min_dis_to_obs(
            [data.x1[1], data.x2[1], data.x3[1], data.x4[1]],
            xo
        )
        data.jl[1] = check_JointLimitation(q₀)
        data.F_distur[1] = zero(u₀)
        data.Fc[1] = zero(u₀)



    end





    # xi = [x1, x2, x3, x4]


    # xs = [[0.0, 0.0]]
    # for x in xi
    #     push!(xs, x(q))
    # end

    # x, y = split_vec_of_arrays(xs)
    # plot(x, y)

    # scatter!([xd[1]], [xd[2]])




end




@time run_simulation()