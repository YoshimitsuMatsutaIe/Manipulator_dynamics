using YAML
using LinearAlgebra
using Dates

#push!(LOAD_PATH, ".")  # includeの代わり．includeダメ絶対


include("utils.jl")
include("rmp.jl")
include("rmp_tree.jl")
include("static_environment.jl")
include("plot_using_Plots.jl")
include("kinematics.jl")
include("lagrange.jl")




# using .RMP
# using .RMPTree
#using .StaticEnvironment
using .Kinematics: q_neutral, q_min, q_max
#using .Utilis
using .Dynamics


### yaml関係 ###

"""動かない障害物をセット

obs_param : yamlを読んで作ったObsParam_の入ったリスト  
"""
function set_obs(obs_param)
    if isnothing(obs_param)
        return nothing
    else
        obs = Vector{State}()
        for param in obs_param
            p = param["data"] |> keytosymbol
            if param["name"] == "point"
                arg = ObsParam_point(;p...)
            elseif param["name"] == "cylinder"
                arg = ObsParam_cylinder(;p...)
            elseif param["name"] == "plane"
                arg = ObsParam_plane(;p...)
            end
            append!(obs, _set_obs(arg))
        end
    end
    return obs
end


"""静止目標をセット"""
function set_goal(goal_param)
    p = goal_param["data"] |> keytosymbol
    if goal_param["name"] == "static"
        arg = GoalParam_point(;p...)
    end
    _set_goal(arg)
end


"""rmp制御器をセット"""
function set_rmp(rmp_param)
    goal = []
    if rmp_param[1]["collision_avoidance"]["name"] == "Original"
        obs = OriginalRMPCollisionAvoidance{Float64}[]
    elseif rmp_param[1]["collision_avoidance"]["name"] == "fromGDS"
        obs = RMPfromGDSCollisionAvoidance{Float64}[]
    end
    jl = []

    for p in rmp_param
        #println(p)

        # 目標到達
        if haskey(p, "goal_attractor") && !isnothing(p["goal_attractor"])
            pa = p["goal_attractor"]
            _p = pa["data"] |> keytosymbol
            if pa["name"] == "Original"
                push!(goal, OriginalRMPAttractor(;_p...))
            elseif pa["name"] == "fromGDS"
                push!(goal, RMPfromGDSAttractor(;_p...))
            end
        end
        
        # 障害物回避
        if haskey(p, "collision_avoidance") && !isnothing(p["collision_avoidance"])
            po = p["collision_avoidance"]
            _p = po["data"] |> keytosymbol
            if po["name"] == "Original"
                push!(obs, OriginalRMPCollisionAvoidance(;_p...))
            elseif po["name"] == "fromGDS"
                push!(obs, RMPfromGDSCollisionAvoidance(;_p...))
            end
        end

        # ジョイント制限回避
        if haskey(p, "joint_limit_avoidance")
            pjl = p["joint_limit_avoidance"]
            _p = pjl["data"] |> keytosymbol
            if pjl["name"] == "Original"
                push!(jl, OriginalJointLimitAvoidance(;_p...))
            elseif pjl["name"] == "fromGDS"
                push!(jl, RMPfromGDSJointLimitAvoidance(;_p...))
            end
        end
    end

    return (
        joint_limit_avoidance=jl[1],
        attractor=goal[1],
        obs_avoidance=obs,
    )
end


"""ジョイント制限を守れてるかチェック"""
function check_JointLimitation(q::Vector{T}) where T
    return q_min .< q .< q_max
end


"""シミュレーションのデータ

t : 時刻  
q : ジョイントベクトルのリスト  
dq :  
ddq :  
desired_ddq :  
u :  





"""
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
    jl::Vector{BitVector}
end



"""
オイラー法で（質量考えずに）シミュレーション
"""
function whithout_mass(
    q₀::Vector{T}, dq₀::Vector{T}, TIME_SPAN::T, Δt::T, obs, goal, rmp_param
) where T

    t = range(0.0, TIME_SPAN, step = Δt)  # 時間軸

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
        Vector{Vector{State{T}}}(undef, length(t)),
        Vector{BitVector}(undef, length(t))
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
    data.jl[1] = check_JointLimitation(q₀)
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
        data.jl[i+1] = check_JointLimitation(data.q[i+1])
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
質量を考えてシミュレーション実行

ルンゲクッタを使用  
"""
function with_mass(
    q₀::Vector{T}, dq₀::Vector{T}, TIME_SPAN::T, Δt::T, obs, goal, rmp_param
) where T

    t = range(0.0, TIME_SPAN, step = Δt)  # 時間軸

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
        Vector{Vector{State{T}}}(undef, length(t)),
        Vector{BitVector}(undef, length(t))
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
    data.jl[1] = check_JointLimitation(q₀)
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
        data.jl[i+1] = check_JointLimitation(data.q[i+1])
    end

    data
end







"""ひとまずシミュレーションやってみｓる"""
function run_simulation(
    isWithMas::Bool, TIME_SPAN::T, Δt::T,
    rmps, obs, goal, t
) where T

    q₀ = q_neutral
    dq₀ = zeros(T, 7)

    if isWithMas
        data = with_mass(q₀, dq₀, TIME_SPAN, Δt, obs, goal, rmps)
    else
        data = whithout_mass(q₀, dq₀, TIME_SPAN, Δt, obs, goal, rmps)
    end
    
    plot_simulation_data(data, t)
    return data
end




"""ランナー

設定を読み込みシミュレーションを実行  
"""
function runner(config, path)
    params = YAML.load_file(config)
    sim_param = params["sim_param"]
    rmp_param = params["rmp_param"]
    env_param = params["env_param"]
    obs = set_obs(env_param["obstacle"])
    goal = set_goal(env_param["goal"])
    rmps = set_rmp(rmp_param)
    
    data= run_simulation(
        sim_param["isWithMass"],
        sim_param["TIME_SPAN"],
        sim_param["TIME_INTERVAL"],
        rmps,
        obs,
        goal,
        path
    )
    return data
end




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



config = "./config/sice.yaml"  # シミュレーション設定のパス
path = get_time_string()  # 実行時のデータ保存パス

println("hoge...")
@time data = runner(config, path)
println("hoge!")

@time make_animation(data, path)

