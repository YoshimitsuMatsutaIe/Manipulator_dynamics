
"""
sice用のモジュール
"""
module SiceKinematics


const l1 = 1.0
const l2 = 1.0
const l3 = 1.0
const l4 = 1.0



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

function J1_dot(q::Vector{T}, dq::Vector{T}) where T
    [
        -dq[1]*l1*cos(q[1]) - l2*(dq[1] + dq[2])*cos(q[1] + q[2]) -l2*(dq[1] + dq[2])*cos(q[1] + q[2]) 0 0
        -dq[1]*l1*sin(q[1]) - l2*(dq[1] + dq[2])*sin(q[1] + q[2]) -l2*(dq[1] + dq[2])*sin(q[1] + q[2]) 0 0
    ]
end

function J2_dot(q::Vector{T}, dq::Vector{T}) where T
    [
        -dq[1]*l1*cos(q[1]) - l2*(dq[1] + dq[2])*cos(q[1] + q[2]) - l3*(dq[1] + dq[2] + dq[3])*cos(q[1] + q[2] + q[3]) -l2*(dq[1] + dq[2])*cos(q[1] + q[2]) - l3*(dq[1] + dq[2] + dq[3])*cos(q[1] + q[2] + q[3]) -l3*(dq[1] + dq[2] + dq[3])*cos(q[1] + q[2] + q[3]) 0
        -dq[1]*l1*sin(q[1]) - l2*(dq[1] + dq[2])*sin(q[1] + q[2]) - l3*(dq[1] + dq[2] + dq[3])*sin(q[1] + q[2] + q[3]) -l2*(dq[1] + dq[2])*sin(q[1] + q[2]) - l3*(dq[1] + dq[2] + dq[3])*sin(q[1] + q[2] + q[3]) -l3*(dq[1] + dq[2] + dq[3])*sin(q[1] + q[2] + q[3]) 0
    ]
end

function J3_dot(q::Vector{T}, dq::Vector{T}) where T
    [
        -dq[1]*l1*cos(q[1]) - l2*(dq[1] + dq[2])*cos(q[1] + q[2]) - l3*(dq[1] + dq[2] + dq[3])*cos(q[1] + q[2] + q[3]) - l4*(dq[1] + dq[2] + dq[3] + dq[4])*cos(q[1] + q[2] + q[3] + q[4]) -l2*(dq[1] + dq[2])*cos(q[1] + q[2]) - l3*(dq[1] + dq[2] + dq[3])*cos(q[1] + q[2] + q[3]) - l4*(dq[1] + dq[2] + dq[3] + dq[4])*cos(q[1] + q[2] + q[3] + q[4]) -l3*(dq[1] + dq[2] + dq[3])*cos(q[1] + q[2] + q[3]) - l4*(dq[1] + dq[2] + dq[3] + dq[4])*cos(q[1] + q[2] + q[3] + q[4]) -l4*(dq[1] + dq[2] + dq[3] + dq[4])*cos(q[1] + q[2] + q[3] + q[4])
        -dq[1]*l1*sin(q[1]) - l2*(dq[1] + dq[2])*sin(q[1] + q[2]) - l3*(dq[1] + dq[2] + dq[3])*sin(q[1] + q[2] + q[3]) - l4*(dq[1] + dq[2] + dq[3] + dq[4])*sin(q[1] + q[2] + q[3] + q[4]) -l2*(dq[1] + dq[2])*sin(q[1] + q[2]) - l3*(dq[1] + dq[2] + dq[3])*sin(q[1] + q[2] + q[3]) - l4*(dq[1] + dq[2] + dq[3] + dq[4])*sin(q[1] + q[2] + q[3] + q[4]) -l3*(dq[1] + dq[2] + dq[3])*sin(q[1] + q[2] + q[3]) - l4*(dq[1] + dq[2] + dq[3] + dq[4])*sin(q[1] + q[2] + q[3] + q[4]) -l4*(dq[1] + dq[2] + dq[3] + dq[4])*sin(q[1] + q[2] + q[3] + q[4])
    ]
end

end