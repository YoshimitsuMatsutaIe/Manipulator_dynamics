using CPUTime
using Plots
using LinearAlgebra
#using Parameters

include("../utils.jl")

include("./rmp.jl")
include("../kinematics/old.jl")


# function jisaku_solve_euler(dx, x₀, t_span, Δt)
#     """オイラー法
#     ・参考にしました: https://twitter.com/genkuroki/status/1301832571131633665/photo/1
#     """

#     t = range(t_span..., step = Δt)  # 時間軸
#     x = Vector{typeof(x₀)}(undef, length(t))  # 解を格納する1次元配列

#     x[1] = x₀  # 初期値
#     for i in 1:length(x)-1
#         x[i+1] = x[i] + dx(t[i], x[i])Δt
#     end

#     t, x
# end


# function jisaku_solve_RungeKutta(dx, x₀, t_span, Δt)
#     """ルンゲクッタ法（4次）"""

#     t = range(t_span..., step = Δt)  # 時間軸
#     x = Vector{typeof(x₀)}(undef, length(t))  # 解を格納する1次元配列

#     x[1] = x₀  # 初期値
#     for i in 1:length(x)-1
#         k₁ = dx(t[i], x[i])
#         k₂ = dx(t[i]+Δt/2, x[i]+k₁*Δt/2)
#         k₃ = dx(t[i]+Δt/2, x[i]+k₂*Δt/2)
#         k₄ = dx(t[i]+Δt, x[i]+k₃*Δt)
#         x[i+1] = x[i] + (k₁ + 2k₂ + 2k₃ +k₄)Δt/6
#     end

#     t, x
# end

"""点の位置と速度"""
mutable struct State{T}
    x::Vector{T}
    dx::Vector{T}
end


"""ノード（今回は制御点+ジョイント位置点）"""
mutable struct Node{T}
    x::Vector{T}  # 位置
    dx::Vector{T}
    Jax::Matrix{T}
    Jay::Matrix{T}
    Jaz::Matrix{T}
    Jo::Matrix{T}
end


"""システムの状態"""
mutable struct SysthemState{T}
    obs::Vector{State{T}}
end


"""加速度指令を計算（resolve演算結果を返す）"""
function calc_ddq(
    nodes::Vector{Vector{Node{T}}},
    goal::State{T},
    obs::Vector{State{T}}
) where T

    root_f = Vector{T}(undef, 7)
    root_M = Matrix{T}(undef, 7, 7)

    attractor = OriginalRMPAttractor(2.0, 10.0, 0.15, 1.0, 1.0, 5.0)
    obs_avovidance = OriginalRMPCollisionAvoidance(0.2, 1.0, 0.5, 0.5, 1.0)
    joint_limit_avoidance = OriginalJointLimitAvoidance(0.05, 0.1, 0.7)

    for i in 1:9
        if i == 1
            _f, _M = get_natural(
                joint_limit_avoidance, nodes[i][1].x, nodes[i][1].dx, q_max, q_min
            )
            root_f += _f
            root_M += _M
        elseif i == 9
            _f, _M = get_natural(
                attractor, nodes[i][1].x, nodes[i][1].dx, goal.x
            )
            _pulled_f, _pulled_M = pullbacked_rmp(_f, _M, nodes[i][1].J,)
            root_f += _pulled_f
            root_M += _pulled_M
        else
            for j in 1:length(nodes[i])
                for k in 1:length(obs)
                    _f, _M = get_natural(
                        obs_avovidance, nodes[i][j].x, nodes[i][j].dx, obs[k].x
                    )
                    _pulled_f, _pulled_M = pullbacked_rmp(_f, _M, nodes[i][j].J,)
                    root_f += _pulled_f
                    root_M += _pulled_M
                end
            end
        end
    end

    ddq = pinv(pulled_M) * pulled_f
    return ddq
end


function update_nodes(nodes::Vector{Vector{Node{T}}}, q::Vector{T}, dq::Vector{T}) where T

end


"""ひとまずシミュレーションやってみｓる"""
function run_simulation()

    TIME_SPAN = 10.0
    Δt = 0.01

    goal = State([1.0 ,1.0, 1.0], [0.0, 0.0, 0.0])
    obs = [
        State([1.0 ,1.0, 1.0], [0.0, 0.0, 0.0])
    ]

    println(typeof(obs))



end



@time run_simulation()