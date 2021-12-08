
"""
sice用のモジュール
"""
module SiceKinematics

# export x1, x2, x3, x4
# export J1, J2, J3, J4
# export J1_dot, J2_dot, J3_dot, J4_dot

export calc_x, calc_J, calc_dJ


const l1 = 1.0
const l2 = 1.0
const l3 = 1.0
const l4 = 1.0

function x1(q)
    [
        l1*cos(q[1])
        l1*sin(q[1])
    ]
end

function x2(q)
    [
        l1*cos(q[1]) + l2*cos(q[1] + q[2])
        l1*sin(q[1]) + l2*sin(q[1] + q[2])
    ]
end

function x3(q)
    [
        l1*cos(q[1]) + l2*cos(q[1] + q[2]) + l3*cos(q[1] + q[2] + q[3])
        l1*sin(q[1]) + l2*sin(q[1] + q[2]) + l3*sin(q[1] + q[2] + q[3])
    ]
end

function x4(q)
    [
        l1*cos(q[1]) + l2*cos(q[1] + q[2]) + l3*cos(q[1] + q[2] + q[3]) + l4*cos(q[1] + q[2] + q[3] + q[4])
        l1*sin(q[1]) + l2*sin(q[1] + q[2]) + l3*sin(q[1] + q[2] + q[3]) + l4*sin(q[1] + q[2] + q[3] + q[4])
    ]
end



function J1(q::Vector{T}) where T
    [
        -l1*sin(q[1]) 0 0 0
        l1*cos(q[1]) 0 0 0
    ]
end

function J2(q::Vector{T}) where T
    [
        -l1*sin(q[1]) - l2*sin(q[1] + q[2]) -l2*sin(q[1] + q[2]) 0 0
        l1*cos(q[1]) + l2*cos(q[1] + q[2]) l2*cos(q[1] + q[2]) 0 0
    ]
end

function J3(q::Vector{T}) where T
    [
        -l1*sin(q[1]) - l2*sin(q[1] + q[2]) - l3*sin(q[1] + q[2] + q[3]) -l2*sin(q[1] + q[2]) - l3*sin(q[1] + q[2] + q[3]) -l3*sin(q[1] + q[2] + q[3]) 0
        l1*cos(q[1]) + l2*cos(q[1] + q[2]) + l3*cos(q[1] + q[2] + q[3]) l2*cos(q[1] + q[2]) + l3*cos(q[1] + q[2] + q[3]) l3*cos(q[1] + q[2] + q[3]) 0
    ]
end

function J4(q::Vector{T}) where T
    [
        -l1*sin(q[1]) - l2*sin(q[1] + q[2]) - l3*sin(q[1] + q[2] + q[3]) - l4*sin(q[1] + q[2] + q[3] + q[4]) -l2*sin(q[1] + q[2]) - l3*sin(q[1] + q[2] + q[3]) - l4*sin(q[1] + q[2] + q[3] + q[4]) -l3*sin(q[1] + q[2] + q[3]) - l4*sin(q[1] + q[2] + q[3] + q[4]) -l4*sin(q[1] + q[2] + q[3] + q[4])
        l1*cos(q[1]) + l2*cos(q[1] + q[2]) + l3*cos(q[1] + q[2] + q[3]) + l4*cos(q[1] + q[2] + q[3] + q[4]) l2*cos(q[1] + q[2]) + l3*cos(q[1] + q[2] + q[3]) + l4*cos(q[1] + q[2] + q[3] + q[4]) l3*cos(q[1] + q[2] + q[3]) + l4*cos(q[1] + q[2] + q[3] + q[4]) l4*cos(q[1] + q[2] + q[3] + q[4])
    ]
end


function J1_dot(q::Vector{T}, dq::Vector{T}) where T
    [
        -dq[1]*l1*cos(q[1]) 0 0 0
        -dq[1]*l1*sin(q[1]) 0 0 0
    ]
end

function J2_dot(q::Vector{T}, dq::Vector{T}) where T
    [
        -dq[1]*l1*cos(q[1]) - l2*(dq[1] + dq[2])*cos(q[1] + q[2]) -l2*(dq[1] + dq[2])*cos(q[1] + q[2]) 0 0
        -dq[1]*l1*sin(q[1]) - l2*(dq[1] + dq[2])*sin(q[1] + q[2]) -l2*(dq[1] + dq[2])*sin(q[1] + q[2]) 0 0
    ]
end

function J3_dot(q::Vector{T}, dq::Vector{T}) where T
    [
        -dq[1]*l1*cos(q[1]) - l2*(dq[1] + dq[2])*cos(q[1] + q[2]) - l3*(dq[1] + dq[2] + dq[3])*cos(q[1] + q[2] + q[3]) -l2*(dq[1] + dq[2])*cos(q[1] + q[2]) - l3*(dq[1] + dq[2] + dq[3])*cos(q[1] + q[2] + q[3]) -l3*(dq[1] + dq[2] + dq[3])*cos(q[1] + q[2] + q[3]) 0
        -dq[1]*l1*sin(q[1]) - l2*(dq[1] + dq[2])*sin(q[1] + q[2]) - l3*(dq[1] + dq[2] + dq[3])*sin(q[1] + q[2] + q[3]) -l2*(dq[1] + dq[2])*sin(q[1] + q[2]) - l3*(dq[1] + dq[2] + dq[3])*sin(q[1] + q[2] + q[3]) -l3*(dq[1] + dq[2] + dq[3])*sin(q[1] + q[2] + q[3]) 0
    ]
end

function J4_dot(q::Vector{T}, dq::Vector{T}) where T
    [
        -dq[1]*l1*cos(q[1]) - l2*(dq[1] + dq[2])*cos(q[1] + q[2]) - l3*(dq[1] + dq[2] + dq[3])*cos(q[1] + q[2] + q[3]) - l4*(dq[1] + dq[2] + dq[3] + dq[4])*cos(q[1] + q[2] + q[3] + q[4]) -l2*(dq[1] + dq[2])*cos(q[1] + q[2]) - l3*(dq[1] + dq[2] + dq[3])*cos(q[1] + q[2] + q[3]) - l4*(dq[1] + dq[2] + dq[3] + dq[4])*cos(q[1] + q[2] + q[3] + q[4]) -l3*(dq[1] + dq[2] + dq[3])*cos(q[1] + q[2] + q[3]) - l4*(dq[1] + dq[2] + dq[3] + dq[4])*cos(q[1] + q[2] + q[3] + q[4]) -l4*(dq[1] + dq[2] + dq[3] + dq[4])*cos(q[1] + q[2] + q[3] + q[4])
        -dq[1]*l1*sin(q[1]) - l2*(dq[1] + dq[2])*sin(q[1] + q[2]) - l3*(dq[1] + dq[2] + dq[3])*sin(q[1] + q[2] + q[3]) - l4*(dq[1] + dq[2] + dq[3] + dq[4])*sin(q[1] + q[2] + q[3] + q[4]) -l2*(dq[1] + dq[2])*sin(q[1] + q[2]) - l3*(dq[1] + dq[2] + dq[3])*sin(q[1] + q[2] + q[3]) - l4*(dq[1] + dq[2] + dq[3] + dq[4])*sin(q[1] + q[2] + q[3] + q[4]) -l3*(dq[1] + dq[2] + dq[3])*sin(q[1] + q[2] + q[3]) - l4*(dq[1] + dq[2] + dq[3] + dq[4])*sin(q[1] + q[2] + q[3] + q[4]) -l4*(dq[1] + dq[2] + dq[3] + dq[4])*sin(q[1] + q[2] + q[3] + q[4])
    ]
end

"""位置を全部計算"""
function calc_x(q::Vector{T}) where T
    x1(q), x2(q), x3(q), x4(q)
end

"""ヤコビ行列を全部計算"""
function calc_J(q::Vector{T}) where T
    return J1(q), J2(q), J3(q), J4(q)
end

"""ヤコビ行列の時間微分を全部計算"""
function calc_dJ(q::Vector{T}, dq::Vector{T}) where T
    J1_dot(q, dq), J2_dot(q, dq), J3_dot(q, dq), J4_dot(q, dq)
end


end