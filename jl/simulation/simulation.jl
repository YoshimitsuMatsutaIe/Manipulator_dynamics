using YAML
using LinearAlgebra
using Dates

include("../utils.jl")
include("../rmp/rmp.jl")
include("../rmp/rmp_tree.jl")
include("static_environment.jl")
include("plot_using_Plots.jl")
include("../kinematics/kinematics.jl")
include("../dynamics/lagrange.jl")

# using .RMP
# using .RMPTree
# using .StaticEnvironment
# using .Kinematics: q_neutral
# using .Utilis
using .Dynamics


"""シミュレーションのデータ"""
struct Data{T}
    t::StepRangeLen{T}  # 時刻
    q::Vector{Vector{T}}  # ジョイントベクトル
    dq::Vector{Vector{T}}  # ジョイント角速度ベクトル
    ddq::Vector{Vector{T}}  # ジョイント角加速度ベクトル
    desired_ddq::Vector{Vector{T}}  # 制御入力ベクトル
    u::Vector{Vector{T}}  # トルクベクトル
    error::Vector{T}  # eeと目標位置との誤差
    dis_to_obs::Vector{T}
    nodes::Vector{Vector{Vector{Node{T}}}}  # ノード
    goal::Vector{State{T}}
    obs::Vector{Vector{State{T}}}
end



"""
オイラー法で（質量考えずに）シミュレーション
"""
function whithout_mass(q₀::Vector{T}, dq₀::Vector{T}, TIME_SPAN::T, Δt::T, obs) where T


    rmp_param = (
        joint_limit_avoidance=OriginalJointLimitAvoidance(0.05, 0.1, 0.7),
        attractor=RMPfromGDSAttractor(10.0, 20.0, 0.15, 2.0, 2.0, 1.0, 0.01, 0.15, 1e-5),
        obs_avoidance=(
            RMPfromGDSCollisionAvoidance(0.5, 1.0, 0.1),
            RMPfromGDSCollisionAvoidance(0.5, 1.0, 0.1),
            RMPfromGDSCollisionAvoidance(0.5, 1.0, 0.1),
            RMPfromGDSCollisionAvoidance(0.5, 1.0, 0.1),
            RMPfromGDSCollisionAvoidance(0.5, 1.0, 0.1),
            RMPfromGDSCollisionAvoidance(0.5, 1.0, 0.1),
            RMPfromGDSCollisionAvoidance(0.5, 1.0, 0.1),
            RMPfromGDSCollisionAvoidance(0.5, 1.0, 0.1),
            RMPfromGDSCollisionAvoidance(0.5, 1.0, 0.1),
        )
    )


    t = range(0.0, TIME_SPAN, step = Δt)  # 時間軸
    #goal = State([0.3, -0.75, 1.0], [0.0, 0.0, 0.0])
    goal = State([0.5, -0.5, 1.3], [0.0, 0.0, 0.0])

    # 初期値
    nodes₀ = update_nodes(nothing, q₀, dq₀)

    data = Data(
        t,
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
        [zeros(T, 7) for i in 1:length(t)],
        Vector{T}(undef, length(t)),
        Vector{T}(undef, length(t)),
        Vector{Vector{Vector{Node{T}}}}(undef, length(t)),
        Vector{State{T}}(undef, length(t)),
        Vector{Vector{State{T}}}(undef, length(t))
    )

    data.q[1] = q₀
    data.dq[1] = dq₀
    data.ddq[1] = zeros(T, 7)
    data.desired_ddq[1] = zeros(T, 7)
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

        data.desired_ddq[i+1], data.dis_to_obs[i+1] = calc_desired_ddq(data.nodes[i], rmp_param, data.goal[i], data.obs[i])
        data.ddq[i+1] = data.desired_ddq[i+1]
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
ルンゲクッタ用

q : 関節角度ベクトル  
dq : 関節角速度ベクトル  
u : 入力トルクベクトル  
F : 外力ベクトル  
"""
function dx(q::Vector{T}, dq::Vector{T}, u::Vector{T}, F::Vector{T}) where T
    ddq = calc_real_ddq(u, F, q, dq)
    return (dq = dq , ddq = ddq)
end


"""
トルク考えてシミュレーション実行

ルンゲクッタを私用
"""
function with_mass(q₀::Vector{T}, dq₀::Vector{T}, TIME_SPAN::T, Δt::T, obs) where T


    rmp_param = (
        joint_limit_avoidance=OriginalJointLimitAvoidance(0.05, 0.1, 0.7),
        attractor=RMPfromGDSAttractor(10.0, 20.0, 0.15, 2.0, 2.0, 1.0, 0.01, 0.15, 1e-5),
        obs_avoidance=(
            RMPfromGDSCollisionAvoidance(0.5, 1.0, 0.1),
            RMPfromGDSCollisionAvoidance(0.5, 1.0, 0.1),
            RMPfromGDSCollisionAvoidance(0.5, 1.0, 0.1),
            RMPfromGDSCollisionAvoidance(0.5, 1.0, 0.1),
            RMPfromGDSCollisionAvoidance(0.5, 1.0, 0.1),
            RMPfromGDSCollisionAvoidance(0.5, 1.0, 0.1),
            RMPfromGDSCollisionAvoidance(0.5, 1.0, 0.1),
            RMPfromGDSCollisionAvoidance(0.5, 1.0, 0.1),
            RMPfromGDSCollisionAvoidance(0.5, 1.0, 0.1),
        )
    )


    t = range(0.0, TIME_SPAN, step = Δt)  # 時間軸
    #goal = State([0.3, -0.75, 1.0], [0.0, 0.0, 0.0])
    goal = State([0.5, -0.5, 1.1], [0.0, 0.0, 0.0])

    # 初期値
    nodes₀ = update_nodes(nothing, q₀, dq₀)

    data = Data(
        t,
        Vector{Vector{T}}(undef, length(t)),
        Vector{Vector{T}}(undef, length(t)),
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
    data.desired_ddq[1] = zeros(T, 7)
    data.u[1] = zeros(T, 7)
    data.error[1] = norm(goal.x - nodes₀[9][1].x)
    data.dis_to_obs[1] = 0.0
    data.nodes[1] = nodes₀
    data.goal[1] = goal
    data.obs[1] = obs
    #println(length(obs))


    F = zeros(T, 7)

    # ぐるぐる回す
    for i in 1:length(data.t)-1
        #println("i = ", i)

        data.nodes[i+1] = update_nodes(data.nodes[i], data.q[i], data.dq[i])
        
        
        data.error[i+1] = norm(data.goal[i].x .- data.nodes[i][9][1].x)

        data.desired_ddq[i+1], data.dis_to_obs[i+1] = calc_desired_ddq(data.nodes[i], rmp_param, data.goal[i], data.obs[i])
        data.u[i+1] = calc_torque(data.q[i], data.dq[i], data.desired_ddq[i+1])

        k1 = dx(data.q[i], data.dq[i], data.u[i+1], F)
        k2 = dx(data.q[i] .+ k1.dq .* Δt/2, data.dq[i] .+ k1.ddq .* Δt/2, data.u[i+1], F)
        k3 = dx(data.q[i] .+ k2.dq .* Δt/2, data.dq[i] .+ k2.ddq .* Δt/2, data.u[i+1], F)
        k4 = dx(data.q[i] .+ k3.dq .* Δt, data.dq[i] .+ k3.ddq .* Δt, data.u[i+1], F)

        data.q[i+1] = data.q[i] .+ (k1.dq .+ 2 .* k2.dq .+ 2 .* k3.dq .+ k4.dq) .* Δt/6
        data.dq[i+1] = data.dq[i] .+ (k1.ddq .+ 2 .* k2.ddq .+ 2 .* k3.ddq .+ k4.ddq) .* Δt/6
        data.ddq[i+1] = k1.ddq .+ 2 .* k2.ddq .+ 2 .* k3.ddq .+ k4.ddq
        #ddq[i+1] = zeros(T,7)
        #println("q = ", q[i])
        #println("dq = ", dq[i])
        #println("ddq = ", ddq[i])

        #data.q[i+1] = data.q[i] .+ data.dq[i]*Δt
        #q[i+1] = zeros(T,7)
        #data.dq[i+1] = data.dq[i] .+ data.ddq[i]*Δt
        data.goal[i+1] = goal
        data.obs[i+1] = obs
    end

    data
end







"""ひとまずシミュレーションやってみｓる"""
function run_simulation(TIME_SPAN::T, Δt::T, obs, t) where T

    q₀ = q_neutral
    dq₀ = zeros(T, 7)

    data = whithout_mass(q₀, dq₀, TIME_SPAN, Δt, obs)
    #data = with_mass(q₀, dq₀, TIME_SPAN, Δt, obs)
    
    plot_simulation_data(data, t)
    return data
end





function runner(name, t)
    params = YAML.load_file(name)
    sim_param = params["sim_param"]
    rmp_param = params["rmp_param"]
    env_param = params["env_param"]
    obs = set_obs(env_param["obstacle"])
    #print(obs)
    data= run_simulation(
        sim_param["TIME_SPAN"], sim_param["TIME_INTERVAL"], obs, t
    )
    return data
end


# using profview
# @profview data, fig = runner("./config/use_RMPfromGDS_test.yaml")

#pass = 
pass = "./config/sice.yaml"


"""データ放送のパス

シミュレーション結果を保存するディレクトリのパスを返す  
ないなら新規作成  
"""
function get_time_string()
    _nd = now()  # 時刻取得
    _Y = _nd |> Dates.year |> string
    _M = _nd |> Dates.month |> string
    _D = _nd |> Dates.day |> string
    _h = _nd |> Dates.hour |> string
    _m = _nd |> Dates.minute |> string
    _s = _nd |> Dates.second |> string
    t = _Y * _M * _D * "-" * _h * _m * _s

    linpath = "../result_of_manipulator_dynamics/" * _Y * _M * _D * "/"
    if !ispath(linpath)
        mkdir(linpath)
    end
    path = linpath * t
    return path
end

path = get_time_string()
println("hoge...")
@time data = runner(pass, path)
println("hoge!")

@time make_animation(data, path)

