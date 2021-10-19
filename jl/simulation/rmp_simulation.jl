using CPUTime
using Plots
using LinearAlgebra
#using Parameters


using PyCall
np = pyimport("numpy")

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
            # for j in 1:length(nodes[i])
            #     for k in 1:length(obs)
            #         _f, _M = get_natural(
            #             obs_avovidance, nodes[i][j].x, nodes[i][j].dx, obs[k].x
            #         )
            #         _pulled_f, _pulled_M = pullbacked_rmp(_f, _M, nodes[i][j].Jo,)
            #         root_f += _pulled_f
            #         root_M += _pulled_M
            #     end
            # end
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





function euler_method(q₀::Vector{T}, dq₀::Vector{T}, TIME_SPAN::T, Δt::T) where T

    t = range(0.0, TIME_SPAN, step = Δt)  # 時間軸
    q = Vector{Vector{T}}(undef, length(t))  # 解を格納する1次元配列
    dq = Vector{Vector{T}}(undef, length(t))
    ddq = Vector{Vector{T}}(undef, length(t))
    error = Vector{T}(undef, length(t))

    goal = State([0.3, -0.75, 1.0], [0.0, 0.0, 0.0])
    obs = [
        State([0.25, -0.4, 1.0], [0.0, 0.0, 0.0])
    ]

    # 初期値
    #println("OK0")
    nodes = update_nodes(nothing, q₀, dq₀)

    error[1] = norm(goal.x - nodes[9][1].x)


    #println("OK")
    q[1] = q₀
    dq[1] = dq₀
    ddq[1] = zeros(T, 7)



    for i in 1:length(q)-1
        #println("i = ", i)
        nodes = update_nodes(nodes, q[i], dq[i])
        #println("OK3")
        error[i+1] = norm(goal.x - nodes[9][1].x)

        ddq[i+1] = calc_ddq(nodes, goal, obs)
        #ddq[i+1] = zeros(T,7)
        #println("q = ", q[i])
        #println("dq = ", dq[i])
        println("ddq = ", ddq[i])

        q[i+1] = q[i] + dq[i]*Δt
        #q[i+1] = zeros(T,7)
        dq[i+1] = dq[i] + ddq[i]*Δt
    end

    t, q, dq, ddq, error
end



"""ひとまずシミュレーションやってみｓる"""
function run_simulation(TIME_SPAN::T, Δt::T) where T

    #q₀ = q_neutral
    q₀ = zeros(T, 7)
    dq₀ = zeros(T, 7)

    t, q, dq, ddq, error = euler_method(q₀, dq₀, TIME_SPAN, Δt)

    q1, q2, q3, q4, q5, q6, q7 = split_vec_of_arrays(q)
    fig_q = plot(t, q1, ylabel="q")
    plot!(fig_q, t, q2)
    plot!(fig_q, t, q3)
    plot!(fig_q, t, q4)
    plot!(fig_q, t, q5)
    plot!(fig_q, t, q6)
    plot!(fig_q, t, q7)

    q1, q2, q3, q4, q5, q6, q7 = split_vec_of_arrays(dq)
    fig_dq = plot(t, q1, ylabel="dq")
    plot!(fig_dq, t, q2)
    plot!(fig_dq, t, q3)
    plot!(fig_dq, t, q4)
    plot!(fig_dq, t, q5)
    plot!(fig_dq, t, q6)
    plot!(fig_dq, t, q7)

    q1, q2, q3, q4, q5, q6, q7 = split_vec_of_arrays(ddq)
    fig_ddq = plot(t, q1, ylabel="ddq")
    plot!(fig_ddq, t, q2)
    plot!(fig_ddq, t, q3)
    plot!(fig_ddq, t, q4)
    plot!(fig_ddq, t, q5)
    plot!(fig_ddq, t, q6)
    plot!(fig_ddq, t, q7)

    fig_error = plot(t, error, ylabel="error")

    plot(
        fig_q, fig_dq, fig_ddq, fig_error, layout=(4,1),
        size=(500,1200)
    )

end



@time run_simulation(1.5, 0.01)
