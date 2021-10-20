using CPUTime
using Plots
using LinearAlgebra
#using Parameters

import YAML


include("../utils.jl")

include("./rmp.jl")
include("../kinematics/old.jl")


"""点の位置と速度"""
mutable struct State{T}
    x::Vector{T}
    dx::Vector{T}
end


"""ノード（今回は制御点+ジョイント位置点）"""
mutable struct Node{T}
    x::Vector{T}  # 位置
    dx::Vector{T}
    # Jax::Matrix{T}
    # Jay::Matrix{T}
    # Jaz::Matrix{T}
    Jo::Matrix{T}
end



"""加速度指令を計算（resolve演算結果を返す）"""
function calc_ddq(
    nodes::Vector{Vector{Node{T}}},
    goal::State{T},
    obs::Vector{State{T}}
) where T

    #root_f = Vector{T}(undef, 7)
    #root_M = Matrix{T}(undef, 7, 7)
    root_f = zeros(T, 7)
    root_M = zeros(T, 7, 7)

    attractor = OriginalRMPAttractor(2.0, 10.0, 0.15, 1.0, 1.0, 5.0)
    obs_avovidance = OriginalRMPCollisionAvoidance(0.2, 1.0, 0.5, 0.5, 1.0)
    joint_limit_avoidance = OriginalJointLimitAvoidance(0.05, 0.1, 0.7)

    for i in 1:9
        if i == 1
            # _f, _M = get_natural(
            #     joint_limit_avoidance, nodes[i][1].x, nodes[i][1].dx, q_max, q_min
            # )
            # root_f += _f
            # root_M += _M
        elseif i == 9
            _f, _M = get_natural(
                attractor, nodes[i][1].x, nodes[i][1].dx, goal.x
            )
            _pulled_f, _pulled_M = pullbacked_rmp(_f, _M, nodes[i][1].Jo,)
            root_f += _pulled_f
            root_M += _pulled_M
        else
            for j in 1:length(nodes[i])
                for k in 1:length(obs)
                    _f, _M = get_natural(
                        obs_avovidance, nodes[i][j].x, nodes[i][j].dx, obs[k].x
                    )
                    _pulled_f, _pulled_M = pullbacked_rmp(_f, _M, nodes[i][j].Jo,)
                    root_f += _pulled_f
                    root_M += _pulled_M
                end
            end
        end
    end

    #root_f += zeros(T, 7)
    #root_M += zeros(T, 7, 7)

    #println("root_f = ", root_f)
    #println("root_M = ", root_M)

    ddq = pinv(root_M) * root_f

    #ddq = np.linalg.pinv(root_M) * root_f

    #ddq = zeros(T, 7)
    return ddq
end


"""nodeを更新"""
function update_nodes(nodes::Vector{Vector{Node{T}}}, q::Vector{T}, dq::Vector{T}) where T

    #println("not nothing")
    # 更新
    DHparams = update_DHparams(q)
    HTMs_local, HTMs_global = calc_HTMs_local_and_global(DHparams)
    Jax_all, Jay_all, Jaz_all, Jo_all = calc_dHTMs(HTMs_local, HTMs_global)
    Jos_joint_all, Jos_cpoint_all = calc_jacobians(Jax_all, Jay_all, Jaz_all, Jo_all)
    cpoints_x_global, cpoints_dx_global = calc_cpoint_x_and_dx_global(
        HTMs_global, Jos_cpoint_all, dq
    )
    #println("www")
    for i in 1:9
        #println("sinu")
        #println(nodes[2])
        for j in 1:length(nodes[i])
            #println("???")
            if i == 1
                nodes[i][j].x = q
                nodes[i][j].dx = dq
            else
                nodes[i][j].x = cpoints_x_global[i-1][j]
                nodes[i][j].dx = cpoints_dx_global[i-1][j]
                nodes[i][j].Jo = Jos_cpoint_all[i-1][j]
            end
        end
    end
    #println("hh")
    return nodes
end


"""nodeを新規作成"""
function update_nodes(nodes::Nothing, q::Vector{T}, dq::Vector{T}) where T

    #println("nothing")
    # 更新
    DHparams = update_DHparams(q)
    HTMs_local, HTMs_global = calc_HTMs_local_and_global(DHparams)
    Jax_all, Jay_all, Jaz_all, Jo_all = calc_dHTMs(HTMs_local, HTMs_global)
    Jos_joint_all, Jos_cpoint_all = calc_jacobians(Jax_all, Jay_all, Jaz_all, Jo_all)
    cpoints_x_global, cpoints_dx_global = calc_cpoint_x_and_dx_global(
        HTMs_global, Jos_cpoint_all, dq
    )
    nodes = Vector{Vector{Node{T}}}(undef, 9)
    for i in 1:9
        #println("i = ", i)
        if i == 1
            _children_node = Vector{Node{T}}(undef, 1)
            _children_node[1] = Node(q, dq, Matrix{T}(I, 7, 7))
            #println("wow")
        else
            n = length(cpoints_local[i-1])
            _children_node = Vector{Node{T}}(undef, n)
            for j in 1:n
                _children_node[j] = Node(
                    cpoints_x_global[i-1][j],
                    cpoints_dx_global[i-1][j],
                    Jos_cpoint_all[i-1][j]
                )
            end
        end
        nodes[i] = _children_node
    end
    return nodes
end


"""シミュレーションのデータ"""
mutable struct Data{T}
    t::StepRangeLen{T}  # 時刻
    q::Vector{Vector{T}}  # ジョイントベクトル
    dq::Vector{Vector{T}}  # ジョイント角速度ベクトル
    ddq::Vector{Vector{T}}  # 制御入力ベクトル
    error::Vector{T}  # eeと目標位置との誤差
    nodes::Vector{Vector{Vector{Node{T}}}}  # ノード
    goal::Vector{State{T}}
    obs::Vector{Vector{State{T}}}
end


"""
オイラー法でシミュレーション
"""
function euler_method(q₀::Vector{T}, dq₀::Vector{T}, TIME_SPAN::T, Δt::T) where T

    t = range(0.0, TIME_SPAN, step = Δt)  # 時間軸
    goal = State([0.3, -0.75, 1.0], [0.0, 0.0, 0.0])
    obs = [
        State([0.25, -0.4, 1.0], [0.0, 0.0, 0.0])
    ]

    # 初期値
    nodes₀ = update_nodes(nothing, q₀, dq₀)

    data = Data(
        t,
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{T}(undef, length(t)),
        Vector{Vector{Vector{Node{T}}}}(undef, length(t)),
        Vector{State{T}}(undef, length(t)),
        Vector{Vector{State{T}}}(undef, length(t))
    )

    data.q[1] = q₀
    data.dq[1] = dq₀
    data.ddq[1] = zeros(T, 7)
    data.error[1] = norm(goal.x - nodes₀[9][1].x)
    data.nodes[1] = nodes₀
    data.goal[1] = goal
    data.obs[1] = obs

    # ぐるぐる回す
    for i in 1:length(data.t)-1
        #println("i = ", i)
        data.nodes[i+1] = update_nodes(data.nodes[i], data.q[i], data.dq[i])
        #println("OK3")
        data.error[i+1] = norm(data.goal[i].x - data.nodes[i][9][1].x)

        data.ddq[i+1] = calc_ddq(data.nodes[i], data.goal[i], data.obs[i])
        #ddq[i+1] = zeros(T,7)
        #println("q = ", q[i])
        #println("dq = ", dq[i])
        #println("ddq = ", ddq[i])

        data.q[i+1] = data.q[i] + data.dq[i]*Δt
        #q[i+1] = zeros(T,7)
        data.dq[i+1] = data.dq[i] + data.ddq[i]*Δt
        data.goal[i+1] = goal
        data.obs[i+1] = obs
    end

    data
end




"""ひとまずシミュレーションやってみｓる"""
function run_simulation(TIME_SPAN::T, Δt::T) where T

    q₀ = q_neutral
    dq₀ = zeros(T, 7)

    data = euler_method(q₀, dq₀, TIME_SPAN, Δt)

    q1, q2, q3, q4, q5, q6, q7 = split_vec_of_arrays(data.q)
    fig_q = plot(data.t, q1, ylabel="q")
    plot!(fig_q, data.t, q2)
    plot!(fig_q, data.t, q3)
    plot!(fig_q, data.t, q4)
    plot!(fig_q, data.t, q5)
    plot!(fig_q, data.t, q6)
    plot!(fig_q, data.t, q7)

    q1, q2, q3, q4, q5, q6, q7 = split_vec_of_arrays(data.dq)
    fig_dq = plot(data.t, q1, ylabel="dq")
    plot!(fig_dq, data.t, q2)
    plot!(fig_dq, data.t, q3)
    plot!(fig_dq, data.t, q4)
    plot!(fig_dq, data.t, q5)
    plot!(fig_dq, data.t, q6)
    plot!(fig_dq, data.t, q7)

    q1, q2, q3, q4, q5, q6, q7 = split_vec_of_arrays(data.ddq)
    fig_ddq = plot(data.t, q1, ylabel="ddq")
    plot!(fig_ddq, data.t, q2)
    plot!(fig_ddq, data.t, q3)
    plot!(fig_ddq, data.t, q4)
    plot!(fig_ddq, data.t, q5)
    plot!(fig_ddq, data.t, q6)
    plot!(fig_ddq, data.t, q7)

    fig_error = plot(data.t, data.error, ylabel="error", ylims=(0.0,))

    fig = plot(
        fig_q, fig_dq, fig_ddq, fig_error, layout=(4,1),
        size=(500,1200)
    )

    return data, fig
end



function make_animation(q, dq, nodes, goal=nothing, obs=nothing)
    anim = Animation()
    @gif for i in 1:10:length(t)
        _fig=draw_arm(q[i], dq[i], goal, obs)
        frame(anim,_fig)
    end
    gif(anim, "test5.gif", fps = 12)
end



@time data, fig = run_simulation(1.0, 0.01)
plot(fig)

# @time t, q, dq, ddq, error, fig, fig2= run_simulation(5.0, 0.01)
# make_animation(q, dq)