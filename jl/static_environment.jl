"""環境いろいろ"""



include("rmp_tree.jl")

#using .RMPTree: State

export ObsParam_point
export ObsParam_cylinder
export set_obs


using Random
using Parameters

include("utils.jl")

sd = Random.seed!(123)

### 目標位置 ###

"""点目標"""
@with_kw struct GoalParam_point{T}
    x::T
    y::T
    z::T
end


function _set_goal(p::GoalParam_point{T}) where T
    State(
        [p.x, p.y, p.z], zeros(T, 3)
    )
end



### 障害物 ###
@with_kw struct ObsParam_field{T, U}
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
@with_kw struct ObsParam_point{T}
    x::T
    y::T
    z::T
end


"""障害物（点）を設置"""
function _set_obs(p::ObsParam_point{T}) where T
    obs = Vector{State{T}}(undef, 1)
    obs[1] = State(
        [p.x, p.y, p.z], zeros(T, 3)
    ) 
    return obs
end

"""
障害物（面）のパラメータ
"""
@with_kw struct ObsParam_plane{T, U}
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
    R = rotate(p.alpha, p.beta, p.gamma)
    t = [p.x, p.y, p.z]

    for i in 1:p.n
        X = [
            (rand(sd, T) - 0.5) * p.lx
            (rand(sd, T) - 0.5) * p.ly
            0.0
        ]

        obs[i] = State(
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
@with_kw struct ObsParam_sphere{T, U}
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
        obs[i] = State(
            X, zeros(T, 3)
        )
    end
    return obs
end


"""
障害物（円筒）のパラメータ

r : 
"""
@with_kw struct ObsParam_cylinder{T, U}
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
    R = rotate(p.alpha, p.beta, p.gamma)
    t = [p.x, p.y, p.z]
    for i in 1:p.n
        theta = 2π * rand(sd, T)
        X = [
            p.r * cos(theta)
            p.r * sin(theta)
            (p.L)*(rand(sd, T) - 1/2)
        ]
        obs[i] = State(
            R * X + t,
            zeros(T, 3)
        )
    end
    return obs
end



