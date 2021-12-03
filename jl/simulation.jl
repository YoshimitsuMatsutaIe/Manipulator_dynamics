using YAML
using LinearAlgebra
using Dates
using Random
using Parameters

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
        obs = Vector{State{Float64}}()
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
    return (q .< q_min) .| (q_max .< q)
end


"""制御点と障害物点の最小距離を計算"""
function calc_min_dis_to_obs(
    nodes::Vector{Vector{Node{T}}}, obs::Vector{State{T}}
    ) where T
    dis_to_obs = 1.0  # 障害物への最小近接距離
    
    _temp_dis_to_obs = Vector{T}(undef, length(obs))
    
    for i in 2:9
        for j in 1:length(nodes[i])
            for k in 1:length(obs)
                _temp_dis_to_obs[k] = norm(obs[k].x .- nodes[i][j].x)  # 障害物との距離
            end
            _d = minimum(_temp_dis_to_obs)
            #println(dis_to_obs, _d)
            if i == 1
                dis_to_obs = _d
            elseif dis_to_obs > _d
                #println("k")
                dis_to_obs = _d
            end
        end
        
    end
    return dis_to_obs
end


"""シミュレーションのデータ

t : 時刻  
q : ジョイントベクトルのリスト  
dq :  
ddq :  
desired_ddq :  
u :  
"""
@with_kw struct Data{T}
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


"""初期に行うやつ"""
function init_Data(
    q₀::Vector{T}, dq₀::Vector{T}, TIME_SPAN::T, Δt::T, obs, goal,
    ;isWithMass::Bool,
    ) where T
    
    t = range(0.0, TIME_SPAN, step = Δt)  # 時間軸

    if isWithMass
        _u = Vector{Vector{T}}(undef, length(t))
    else
        _u = [zeros(T, 7) for i in 1:length(t)]
    end
    
    data = Data(
        t = t,
        q = Vector{Vector{T}}(undef, length(t)),
        dq = Vector{Vector{T}}(undef, length(t)),
        ddq = Vector{Vector{T}}(undef, length(t)),
        desired_ddq = Vector{Vector{T}}(undef, length(t)),
        u = _u,
        error = Vector{T}(undef, length(t)),
        dis_to_obs = Vector{T}(undef, length(t)),
        nodes = Vector{Vector{Vector{Node{T}}}}(undef, length(t)),
        goal = Vector{State{T}}(undef, length(t)),
        obs = Vector{Vector{State{T}}}(undef, length(t)),
        jl = Vector{BitVector}(undef, length(t)),
    )

    # 初期値代入
    nodes₀ = update_nodes(nothing, q₀, dq₀)  # 初期ノード

    data.q[1] = q₀
    data.dq[1] = dq₀
    data.ddq[1] = zeros(T, 7)
    data.desired_ddq[1] = zeros(T, 7)
    data.error[1] = norm(goal.x - nodes₀[end][1].x)
    data.dis_to_obs[1] = 0.0
    data.nodes[1] = nodes₀
    data.goal[1] = goal
    data.obs[1] = obs
    data.jl[1] = check_JointLimitation(q₀)

    if isWithMass
        data.u[1] = zeros(T, 7)
    end

    return data
end




"""
オイラー法で（質量考えずに）シミュレーション
"""
function whithout_mass(
    q₀::Vector{T}, dq₀::Vector{T}, TIME_SPAN::T, Δt::T, obs, goal, rmp_param
    ) where T

    data = init_Data(
        q₀, dq₀, TIME_SPAN, Δt, obs, goal, isWithMass=false
        )

    # ぐるぐる回す
    for i in 1:length(data.t)-1
        #println("i = ", i)
        data.nodes[i+1] = update_nodes!(
            nodes = deepcopy(data.nodes[i]),
            q = data.q[i],
            dq = data.dq[i],
            rmp_param = rmp_param,
            goal = data.goal[i],
            obs = data.obs[i],
            )
        #println("OK3")
        data.error[i+1] = norm(data.goal[i].x .- data.nodes[i][end][1].x)
        data.dis_to_obs[i+1] = calc_min_dis_to_obs(data.nodes[i+1], obs)

        data.desired_ddq[i+1] = calc_desired_ddq(
            data.nodes[i]
            )
        data.ddq[i+1] = data.desired_ddq[i+1]

        data.q[i+1] = data.q[i] .+ data.dq[i]*Δt
        #q[i+1] = zeros(T,7)
        data.dq[i+1] = data.dq[i] .+ data.ddq[i]*Δt
        data.goal[i+1] = goal
        data.obs[i+1] = obs
        data.jl[i+1] = check_JointLimitation(data.q[i+1])


        #println(data.nodes[i][1][1].x - data.nodes[i+1][1][1].x)
        #println(data.nodes[i][1] == data.nodes[i+1][1])


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
function dx(;q::Vector{T}, dq::Vector{T}, u::Vector{T}, F::Vector{T}) where T
    ddq = calc_real_ddq(u, F, q, dq)
    return (dq = dq , ddq = ddq)
end


"""
質量を考えてシミュレーション実行

ルンゲクッタを使用  
"""
function with_mass(
    q₀::Vector{T}, dq₀::Vector{T}, TIME_SPAN::T, Δt::T,
    obs::Vector{State{T}}, goal::State{T}, rmp_param
    ) where T

    data = init_Data(
        q₀, dq₀, TIME_SPAN, Δt, obs, goal, isWithMass=true
    )

    ## 一定外乱
    F = zeros(T, 7)  # 外力なし
    #F = rand(T, 7) * 0.001

    # ぐるぐる回す
    for i in 1:length(data.t)-1
        #F = rand(T, 7) * 0.0001
        

        data.nodes[i+1] = update_nodes!(
            nodes = deepcopy(data.nodes[i]),
            q = data.q[i],
            dq = data.dq[i],
            rmp_param = rmp_param,
            goal = data.goal[i],
            obs = data.obs[i],
        )
        data.error[i+1] = norm(data.goal[i].x .- data.nodes[i][end][1].x)
        data.dis_to_obs[i+1] = calc_min_dis_to_obs(data.nodes[i], data.obs[i])

        data.desired_ddq[i+1] = calc_desired_ddq(data.nodes[i])
        data.u[i+1] = calc_torque(data.q[i], data.dq[i], data.desired_ddq[i+1])

        # ルンゲクッタで次の値を決定
        k1 = dx(
            q = data.q[i],
            dq = data.dq[i],
            u = data.u[i+1],
            F = F
            )
        k2 = dx(
            q = data.q[i] .+ k1.dq .* Δt/2,
            dq = data.dq[i] .+ k1.ddq .* Δt/2,
            u = data.u[i+1],
            F = F
            )
        k3 = dx(
            q = data.q[i] .+ k2.dq .* Δt/2,
            dq = data.dq[i] .+ k2.ddq .* Δt/2,
            u = data.u[i+1],
            F = F
            )
        k4 = dx(
            q = data.q[i] .+ k3.dq .* Δt,
            dq = data.dq[i] .+ k3.ddq .* Δt,
            u = data.u[i+1],
            F = F
            )

        data.ddq[i+1] = k1.ddq .+ 2 .* k2.ddq .+ 2 .* k3.ddq .+ k4.ddq
        data.q[i+1] = data.q[i] .+ (k1.dq .+ 2 .* k2.dq .+ 2 .* k3.dq .+ k4.dq) .* Δt/6
        data.dq[i+1] = data.dq[i] .+ data.ddq[i+1] .* Δt/6
        
        
        
        data.goal[i+1] = goal
        data.obs[i+1] = obs
        data.jl[i+1] = check_JointLimitation(data.q[i+1])
    end

    data
end







"""ひとまずシミュレーションやってみｓる"""
function run_simulation(;
    saveRmpData::Bool, isWithMas::Bool, TIME_SPAN::T, Δt::T,
    rmps, obs, goal, save_path
) where T

    q₀ = q_neutral
    dq₀ = zeros(T, 7)

    if isWithMas
        @time data = with_mass(q₀, dq₀, TIME_SPAN, Δt, obs, goal, rmps)
    else
        @time data = whithout_mass(q₀, dq₀, TIME_SPAN, Δt, obs, goal, rmps)
    end
    
    
    @time plot_simulation_data(data, save_path)
    @time plot_rmp(data, save_path)
    @time make_animation(data, save_path)

    #println(typeof(length(data.q)))

    return data
end




"""ランナー

設定を読み込みシミュレーションを実行  
config : yamlファイルのパス  
path : 結果保存先のパス  
"""
function runner(config, save_path)
    params = YAML.load_file(config)
    sim_param = params["sim_param"]
    rmp_param = params["rmp_param"]
    env_param = params["env_param"]
    obs = set_obs(env_param["obstacle"])
    goal = set_goal(env_param["goal"])
    rmps = set_rmp(rmp_param)
    
    data = run_simulation(
        saveRmpData = sim_param["saveRmpData"],
        isWithMas = sim_param["isWithMass"],
        TIME_SPAN = sim_param["TIME_SPAN"],
        Δt = sim_param["TIME_INTERVAL"],
        rmps = rmps,
        obs = obs,
        goal = goal,
        save_path = save_path
    )
    return data
end




"""データ保存先のパス

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
    t = _h * "_" * _m * "_" * _s

    linpath = "../result_of_manipulator_dynamics/" *
    _Y * "_" * _M * "_" * _D * "/" *
    t * "/"

    mkpath(linpath)
    path = linpath
    return path
end



#config = "./config/sice.yaml"  # シミュレーション設定のパス
config = "./config/sice_2.yaml"



path = get_time_string()  # 実行時のデータ保存パス

println("hoge...")
data = runner(config, path)
println("hoge!")

