using CPUTime
using LinearAlgebra
#using Parameters

using YAML


include("../utils.jl")

include("./rmp.jl")
include("../kinematics/old.jl")
include("environment.jl")
include("plot_using_Plots.jl")




"""点の位置と速度"""
mutable struct State{T}
    x::Vector{T}
    dx::Vector{T}
end



function get_x_from_State(obs)
    x = Vector{Vector{typeof(obs[1].x[1])}}(undef, length(obs))
    for i in 1:length(obs)
        x[i] = obs[i].x
    end
    return x
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



#attractor = OriginalRMPAttractor(2.0, 10.0, 0.15, 1.0, 1.0, 5.0)
#const obs_avoidance = OriginalRMPCollisionAvoidance(0.2, 1.0, 0.5, 0.5, 1.0)
const joint_limit_avoidance = OriginalJointLimitAvoidance(0.05, 0.1, 0.7)

const attractor = RMPfromGDSAttractor(5.0, 20.0, 0.15, 2.0, 2.0, 1.0, 0.01, 0.15, 1e-5)
const obs_avoidance = RMPfromGDSCollisionAvoidance(0.1, 1.0, 0.01)



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

    dis_to_obs = 0.0

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
            @. root_f += _pulled_f
            @. root_M += _pulled_M


        else
            for j in 1:length(nodes[i])
                _temp_dis_to_obs = Vector{T}(undef, length(obs))
                for k in 1:length(obs)
                    _temp_dis_to_obs[k] = norm(obs[k].x .- nodes[i][j].x)
                    _f, _M = get_natural(
                        obs_avoidance, nodes[i][j].x, nodes[i][j].dx, obs[k].x
                    )
                    _pulled_f, _pulled_M = pullbacked_rmp(_f, _M, nodes[i][j].Jo,)
                    @. root_f += _pulled_f
                    @. root_M += _pulled_M
                end
                _d = minimum(_temp_dis_to_obs)
                if dis_to_obs < _d
                    dis_to_obs = _d
                end
            end
        end
    end

    #root_f += zeros(T, 7)
    #root_M += zeros(T, 7, 7)

    #println("root_f = ", root_f)
    #println("root_M = ", root_M)


    ddq = pinv(root_M) * root_f
    #println("ddq = ", ddq)
    #ddq = np.linalg.pinv(root_M) * root_f

    #ddq = zeros(T, 7)
    return ddq, dis_to_obs
end


"""nodeを更新"""
function update_nodes(nodes::Vector{Vector{Node{T}}}, q::Vector{T}, dq::Vector{T}) where T

    #println("not nothing")
    # 更新
    HTMs_local, HTMs_global = update_DHparams(q) |> calc_HTMs_local_and_global
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
    HTMs_local, HTMs_global = update_DHparams(q) |> calc_HTMs_local_and_global
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
    dis_to_obs::Vector{T}
    nodes::Vector{Vector{Vector{Node{T}}}}  # ノード
    goal::Vector{State{T}}
    obs::Vector{Vector{State{T}}}
end




"""
オイラー法でシミュレーション
"""
function euler_method(q₀::Vector{T}, dq₀::Vector{T}, TIME_SPAN::T, Δt::T, obs) where T

    t = range(0.0, TIME_SPAN, step = Δt)  # 時間軸
    goal = State([0.3, -0.75, 1.0], [0.0, 0.0, 0.0])

    # 初期値
    nodes₀ = update_nodes(nothing, q₀, dq₀)

    data = Data(
        t,
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{T}(undef, length(t)),
        Vector{T}(undef, length(t)),
        Vector{Vector{Vector{Node{T}}}}(undef, length(t)),
        Vector{State{T}}(undef, length(t)),
        Vector{Vector{State{T}}}(undef, length(t))
    )

    data.q[1] = q₀
    data.dq[1] = dq₀
    data.ddq[1] = zeros(T, 7)
    data.error[1] = norm(goal.x - nodes₀[9][1].x)
    data.dis_to_obs[1] = 0.0
    data.nodes[1] = nodes₀
    data.goal[1] = goal
    data.obs[1] = obs
    #println(length(obs))

    # ぐるぐる回す
    for i in 1:length(data.t)-1
        #println("i = ", i)
        data.nodes[i+1] = update_nodes(data.nodes[i], data.q[i], data.dq[i])
        #println("OK3")
        data.error[i+1] = norm(data.goal[i].x .- data.nodes[i][9][1].x)

        data.ddq[i+1], data.dis_to_obs[i+1] = calc_ddq(data.nodes[i], data.goal[i], data.obs[i])
        #ddq[i+1] = zeros(T,7)
        #println("q = ", q[i])
        #println("dq = ", dq[i])
        #println("ddq = ", ddq[i])

        data.q[i+1] = data.q[i] .+ data.dq[i]*Δt
        #q[i+1] = zeros(T,7)
        data.dq[i+1] = data.dq[i] .+ data.ddq[i]*Δt
        data.goal[i+1] = goal
        data.obs[i+1] = obs
    end

    data
end




"""
ルンゲクッタ．たぶんすごく遅い
"""
function runge_kutta_method(q₀::Vector{T}, dq₀::Vector{T}, TIME_SPAN::T, Δt::T, obs) where T

    t = range(0.0, TIME_SPAN, step = Δt)  # 時間軸
    goal = State([0.3, -0.75, 1.0], [0.0, 0.0, 0.0])

    # 初期値
    nodes₀ = update_nodes(nothing, q₀, dq₀)

    data = Data(
        t,
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{T}(undef, length(t)),
        Vector{T}(undef, length(t)),
        Vector{Vector{Vector{Node{T}}}}(undef, length(t)),
        Vector{State{T}}(undef, length(t)),
        Vector{Vector{State{T}}}(undef, length(t))
    )

    data.q[1] = q₀
    data.dq[1] = dq₀
    data.ddq[1] = zeros(T, 7)
    data.error[1] = norm(goal.x - nodes₀[9][1].x)
    data.dis_to_obs[1] = 0.0
    data.nodes[1] = nodes₀
    data.goal[1] = goal
    data.obs[1] = obs
    #println(length(obs))

    # ぐるぐる回す
    for i in 1:length(data.t)-1
        #println("i = ", i)
        data.nodes[i+1] = update_nodes(data.nodes[i], data.q[i], data.dq[i])
        data.error[i+1] = norm(data.goal[i].x .- data.nodes[i][9][1].x)

        k1 = 

        data.ddq[i+1], data.dis_to_obs[i+1] = calc_ddq(data.nodes[i], data.goal[i], data.obs[i])
        data.q[i+1] = data.q[i] .+ data.dq[i]*Δt
        data.dq[i+1] = data.dq[i] .+ data.ddq[i]*Δt
        data.goal[i+1] = goal
        data.obs[i+1] = obs
    end

    data
end




"""ひとまずシミュレーションやってみｓる"""
function run_simulation(TIME_SPAN::T, Δt::T, obs) where T

    q₀ = q_neutral
    dq₀ = zeros(T, 7)

    data = euler_method(q₀, dq₀, TIME_SPAN, Δt, obs)
    fig = plot_simulation_data(data)
    return data, fig
end





function runner(name)
    params = YAML.load_file(name)
    sim_param = params["sim_param"]
    rmp_param = params["rmp_param"]
    env_param = params["env_param"]
    obs = set_obs(env_param["obstacle"])
    data, fig = run_simulation(
        sim_param["TIME_SPAN"], sim_param["TIME_INTERVAL"], obs
    )
    data, fig
end



@time data, fig = runner("./config/use_RMPfromGDS_test.yaml")
fig
@time make_animation(data)


# @time t, q, dq, ddq, error, fig, fig2= run_simulation(5.0, 0.01)
# make_animation(q, dq)