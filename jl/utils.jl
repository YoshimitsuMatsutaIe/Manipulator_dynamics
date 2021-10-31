"""色々つかうやつ"""

# """
# 様々な便利関数&struct
# """
# module Utilis

# export split_vec_of_arrays
# export rotate_3d
# export State
# export get_x_from_State
# export Node


using StaticArrays
using ArraysOfArrays


function split_vec_of_arrays(u)
    vec.(u) |>
    x ->
    VectorOfSimilarVectors(x).data |>
    transpose |>
    VectorOfSimilarVectors
end

"""3次元回転行列"""
function rotate_3d(a::T, b::T, c::T) where T
    a = deg2rad(a)
    b = deg2rad(b)
    c = deg2rad(c)
    [
        1.0 0.0 0.0
        0.0 cos(a) -sin(a)
        0.0 sin(a) cos(a)
    ] * [
        cos(b) 0.0 sin(b)
        0.0 1.0 0.0
        -sin(b) 0.0 cos(b)
    ] * [
        cos(c) -sin(c) 0.0
        sin(c) cos(c) 0.0
        0.0 0.0 1.0
    ]
end



"""点の位置と速度

x : 位置ベクトル  
dx : 速度ベクトル  
"""
mutable struct State{T}
    x::Vector{T}
    dx::Vector{T}
end


"""
Stateからxだけ取る無駄関数
"""
function get_x_from_State(obs)
    x = Vector{Vector{typeof(obs[1].x[1])}}(undef, length(obs))
    for i in 1:length(obs)
        x[i] = obs[i].x
    end
    return x
end


"""ノード

x : 位置ベクトル  
dx : 速度ベクトル  
Jax : x軸回りの回転軸ベクトルの関節角度ベクトルによるヤコビ行列  
Jay : y軸回りの回転軸ベクトルの関節角度ベクトルによるヤコビ行列  
Jaz : z軸回りの回転軸ベクトルの関節角度ベクトルによるヤコビ行列  
Jo : 位置ベクトルの関節角度ベクトルによるヤコビ行列  
"""
mutable struct Node{T}
    x::Vector{T}  # 位置
    dx::Vector{T}
    Jax::Matrix{T}  # 角度を制御する場合必要
    Jay::Matrix{T}
    Jaz::Matrix{T}
    Jo::Matrix{T}
end


#end