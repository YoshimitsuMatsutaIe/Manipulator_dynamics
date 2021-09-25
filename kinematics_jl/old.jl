using CPUTime
using Plots
using LinearAlgebra

struct DHparam
    alpha::Float64
    a::Float64
    d::Float64
    theta::Float64
end

L = 278e-3
h = 64e-3
H = 1104e-3
L0 = 270.35e-3
L1 = 69e-3
L2 = 364.35e-3
L3 = 69e-3
L4 = 374.29e-3
L5 = 10e-3
L6 = 368.3e-3

q_neutral = [0; -31; 0; 43; 0; 72; 0] * pi/180

q = q_neutral

DHparams = [
    DHparam(0, 0, 0, q[1])
    DHparam(-pi/2, L1, 0, q[2]+pi/2)
    DHparam(pi/2, 0, L2, q[3])
    DHparam(-pi/2, L3, 0, q[4])
    DHparam(pi/2, 0, L4, q[5])
    DHparam(-pi/2, L5, 0, q[6])
    DHparam(pi/2, 0, 0, q[7])
]


T_BR_Wo = [
    -sqrt(2)/2 sqrt(2)/2 0 -L
    -sqrt(2)/2 -sqrt(2)/2 0 -H
    0 0 1 H
    0 0 0 1
]


function update_DHparams(q, DHparams)
    q[2] = q[2] + pi/2
    for i in 1:7
        DHparams[i].theta = q[i]
    end
    DHparams
end




function T(p::DHparam)
    """(i-1)T(i)"""
    [
        cos(p.theta) -sin(p.theta) 0 p.a
        sin(p.theta)*cos(p.alpha) cos(p.theta)*cos(p.alpha) -sin(p.alpha) -p.d*sin(p.alpha)
        sin(p.theta)*sin(p.alpha) cos(p.theta)*sin(p.alpha) cos(p.alpha) p.d*cos(p.alpha)
        0 0 0 1
    ]
end

println(T(DHparams[2]))