"""環境いろいろ"""

using Random

include("../utils.jl")





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
    theta::T
    phi::T
    zeta::T
    n::U
end




function _set_obs(p::ObsParam_point{T}) where T
    obs = Vector{Vector{T}}(undef, 1)
    obs[1] = [p.x, p.y, p.z]
    return obs
end



function _set_obs(p::ObsParam_sphere{T, U}) where {T, U}
    obs = Vector{Vector{T}}(undef, p.n)
    for i in 1:p.n
        theta = 2*rand(T)-1 |> acos
        phi = 2π * rand(T)
        obs[i] = [
            p.r * sin(theta) * cos(phi)
            p.r * sin(theta) * sin(phi)
            p.r * cos(theta)
        ]
    end
    return obs
end

function _set_obs(p::ObsParam_cylinder{T, U}) where {T, U}
