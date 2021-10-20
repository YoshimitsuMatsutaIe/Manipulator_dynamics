"""環境いろいろ"""

using Random

include("../utils.jl")
#include("rmp_simulation.jl")




struct ObsParam_point{T}
    x::T
    y::T
    z::T
end


struct ObsParam_sphere{T, U}
    x::T
    y::T
    z::T
    r::T
    n::U
end


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



function _set_obs(p::ObsParam_point{T}) where T
    obs = Vector{State{T}}(undef, 1)
    obs[1] = State(
        [p.x, p.y, p.z], zeros(T, 3)
    ) 
    return obs
end



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
        obs[i] = State(
            X, zeros(T, 3)
        )
    end
    return obs
end

function _set_obs(p::ObsParam_cylinder{T}) where {T}
    obs = Vector{State{T}}(undef, p.n)
    R = rotate_3d(p.alpha, p.beta, p.gamma)
    t = [p.x, p.y, p.z]
    for i in 1:p.n
        theta = 2π * rand(T)
        X = [
            p.r * cos(theta)
            p.r * sin(theta)
            (p.L)*(rand(T) - 1/2)
        ]
        obs[i] = State(
            R * X + t,
            zeros(T, 3)
        )
    end
    return obs
end


function _set_obs(p::ObsParam_field{T, U}) where {T, U}
    obs = Vector{State{T}}(undef, p.n)
    R = rotate_3d(p.alpha, p.beta, p.gamma)
    t = [p.x, p.y, p.z]
    for i in 1:p.n
        X = [
            p.lx * (rand(T) - 1/2)
            p.ly * (rand(T) - 1/2)
            0.0
        ]
        obs[i] = State(
            R * X + t,
            zeros(T, 3)
        )
    end
    return obs
end



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
            end
            #println(arg)
            append!(obs, _set_obs(arg))
        end
    end
    return obs
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
