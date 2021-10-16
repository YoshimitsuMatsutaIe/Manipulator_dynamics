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
            _f, _M = get_natural(
                joint_limit_avoidance, nodes[i][1].x, nodes[i][1].dx, q_max, q_min
            )
            root_f += _f
            root_M += _M
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

    #println(root_M)
    ddq = pinv(root_M) * root_f
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


"""ひとまずシミュレーションやってみｓる"""
function run_simulation(TIME_SPAN::T, Δt::T) where T



    q₀ = q_neutral
    dq₀ = zeros(T, 7)


    function euler_method(q₀::Vector{T}, dq₀::Vector{T}, TIME_SPAN::T, Δt::T) where T

        t = range(0.0, TIME_SPAN, step = Δt)  # 時間軸
        q = Vector{Vector{T}}(undef, length(t))  # 解を格納する1次元配列
        dq = Vector{Vector{T}}(undef, length(t))  # 解を格納する1次元配列
        
        # 初期値
        #println("OK0")
        nodes = update_nodes(nothing, q₀, dq₀)
        #println("OK")
        q[1] = q₀
        dq[1] = dq₀
        goal = State([0.3, -0.75, 1.0], [0.0, 0.0, 0.0])
        obs = [
            State([0.25, -0.4, 1.0], [0.0, 0.0, 0.0])
        ]

        for i in 1:length(q)-1
            println("i = ", i)
            nodes = update_nodes(nodes, q[i], dq[i])
            #println("OK3")
            ddq = calc_ddq(nodes, goal, obs)
            println(ddq)
            q[i+1] = q[i] + dq[i]*Δt
            dq[i+1] = dq[i] + ddq*Δt
        end
    
        t, q, dq
    end

    t, q, dq = euler_method(q₀, dq₀, TIME_SPAN, Δt)

end



@time t, q, dq = run_simulation(5.0, 0.01)