"""環境いろいろ"""



"""
動かない障害物
"""
module StaticEnvironment

include("../rmp/rmp_tree.jl")

#using .RMPTree: State

export ObsParam_point
export ObsParam_cylinder
export set_obs


using Random

include("../utils.jl")
using .Utilis






struct ObsParam_field{T, U}
    x::T
    y::T
    z::T
    lx::T
    ly::T
    alpha::T
    beta::T
    gamma::T
    n::U
end


"""
障害物（点）のパラメータ
"""
struct ObsParam_point{T}
    x::T
    y::T
    z::T
end


"""障害物（点）を設置"""
function _set_obs(p::ObsParam_point{T}) where T
    obs = Vector{State{T}}(undef, 1)
    obs[1] = Utilis.State(
        [p.x, p.y, p.z], zeros(T, 3)
    ) 
    return obs
end

"""
障害物（面）のパラメータ
"""
struct ObsParam_plane{T, U}
    x::T
    y::T
    z::T
    lx::T
    ly::T
    alpha::T
    beta::T
    gamma::T
    n::U
end


"""
障害物（面）を設置
"""
function _set_obs(p::ObsParam_plane{T, U}) where {T, U}
    obs = Vector{State{T}}(undef, p.n)
    R = Utilis.rotate_3d(p.alpha, p.beta, p.gamma)
    t = [p.x, p.y, p.z]

    for i in 1:p.n
        X = [
            (rand(T) - 0.5) * p.lx
            (rand(T) - 0.5) * p.ly
            0.0
        ]

        obs[i] = Utilis.State(
            R * X + t,
            zeros(T, 3)
        )
    end
    return obs
end


"""
障害物（球）のパラメータ

r : 半径  
n : 障害物点の数  
"""
struct ObsParam_sphere{T, U}
    x::T
    y::T
    z::T
    r::T
    n::U
end


"""
障害物（球）を設置
"""
function _set_obs(p::ObsParam_sphere{T, U}) where {T, U}
    obs = Vector{State{T}}(undef, p.n)
    for i in 1:p.n
        theta = 2*rand(T)-1 |> acos
        phi = 2π * rand(T)
        X = [
            p.r * sin(theta) * cos(phi) + p.x
            p.r * sin(theta) * sin(phi) + p.y
            p.r * cos(theta) + p.z
        ]
        obs[i] = Utilis.State(
            X, zeros(T, 3)
        )
    end
    return obs
end


"""
障害物（円筒）のパラメータ

r : 
"""
struct ObsParam_cylinder{T, U}
    x::T
    y::T
    z::T
    r::T
    L::T
    alpha::T
    beta::T
    gamma::T
    n::U
end


"""障害物（円筒）を設置"""
function _set_obs(p::ObsParam_cylinder{T}) where {T}
    obs = Vector{State{T}}(undef, p.n)
    R = Utilis.rotate_3d(p.alpha, p.beta, p.gamma)
    t = [p.x, p.y, p.z]
    for i in 1:p.n
        theta = 2π * rand(T)
        X = [
            p.r * cos(theta)
            p.r * sin(theta)
            (p.L)*(rand(T) - 1/2)
        ]
        obs[i] = Utilis.State(
            R * X + t,
            zeros(T, 3)
        )
    end
    return obs
end




"""yamelを読んで質量無限大の障害物設置"""
function set_obs(obs_param)
    if isnothing(obs_param)
        return nothing
    else
        obs = Vector{State}()
        for param in obs_param
            p = param["data"]
            if param["name"] == "point"
                arg = ObsParam_point(p["x"], p["y"], p["z"])
            elseif param["name"] == "cylinder"
                arg = ObsParam_cylinder(
                    p["x"], p["y"], p["z"], p["r"], p["L"], p["alpha"], p["beta"], p["gamma"], p["n"]
                )
            elseif param["name"] == "plane"
                arg = ObsParam_plane(
                    p["x"], p["y"], p["z"], p["lx"], p["ly"], p["alpha"], p["beta"], p["gamma"], p["n"]
                )
            end
            #println(arg)
            append!(obs, _set_obs(arg))
        end
    end
    return obs
end




end


# using YAML
# data = YAML.load_file("./config./use_RMPfromGDS_test.yaml")
# obs_param = data["env_param"]["obstacle"]
# println(obs_param)
# o = set_obs(obs_param)
# X = get_x_from_State(o)
# x, y, z = split_vec_of_arrays(X)
# using Plots
# fig = scatter(x, y, z)

using .StaticEnvironment
println(StaticEnvironment.Utilis.rotate_3d(0,0,0))