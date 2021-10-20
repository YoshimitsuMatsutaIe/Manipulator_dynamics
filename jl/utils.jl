"""色々つかうやつ"""

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
function rotate(a::T, b::T, c::T) where T
    [
        1 0 0
        0 cos(a) -sin(a)
        0 sin(a) cos(a)
    ] * [
        cos(b) 0 sin(b)
        0 1 0
        -sin(b) 0 cos(b)
    ] * [
        cos(c) -sin(c) 0
        sin(c) cos(c) 0
        0 0 1
    ]
end


