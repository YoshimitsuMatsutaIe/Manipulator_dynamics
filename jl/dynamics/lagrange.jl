"""ラグランジュ法による入力計算（テスト）
"""

#using CPUTime


"""
ロボットの動力学

運動学とDHパラメータの取り方が異なることに注意  
"""
#module Dynamics


using LinearAlgebra
const n = 7

const Ixx = (
    0.0470910226,
    0.027885975,
    0.0266173355,
    0.0131822787,
    0.0166774282,
    0.0070053791,
    0.0008162135,
)

const Iyy = (
    0.035959884,
    0.020787492,
    0.012480083,
    0.009268520,
    0.003746311,
    0.005527552,
    0.0008735012,
)

const Izz = (
    0.0376697645,
    0.0117520941,
    0.0284435520,
    0.0071158268,
    0.0167545726,
    0.0038760715,
    0.0005494148,
)

const Ixy = (
    -0.0061487003,
    -0.0001882199,
    -0.0039218988,
    -0.0001966341,
    -0.0001865762,
    0.0001534806,
    0.000128440,
)

const Iyz = (
    -0.0007808689,
    0.0020767576,
    -0.001083893,
    0.000745949,
    0.0006473235,
    -0.0002111503,
    0.0001057726,
)

const Ixz = (
    0.0001278755,
    -0.00030096397,
    0.0002927063,
    0.0003603617,
    0.0001840370,
    -0.0004438478,
    0.00018969891,
)

const x_bar = (
    -0.05117,
    0.00269,
    -0.07176,
    0.00159,
    -0.01168,
    0.00697,
    0.005137,
)

const y_bar = (
    0.07908,
    -0.00529,
    0.08149,
    -0.01117,
    0.13111,
    0.006,
    0.0009572,
)

const z_bar = (
    0.00086,
    0.06845,
    0.00132,
    0.02618,
    0.0046,
    0.06048,
    -0.06682,
)

# const r_bars = [
#     -0.05117 0.00269  -0.07176 0.00159  -0.01168 0.00697 0.005137
#     0.07908  -0.00529 0.08149  -0.01117 0.13111  0.006   0.0009572
#     0.00086  0.06845  0.00132  0.02618  0.0046   0.06048 -0.06682
#     0.0      0.0      0.0      0.0      0.0      0.0     0.0
# ]

const m = (
    5.70044,
    3.22698,
    4.31272,
    2.07206,
    2.24665,
    1.60979,
    0.54218,
)

const d = (
    0.2703,
    0.0,
    0.3644,
    0.0,
    0.3743,
    0.0,
    0.2295,
)

const a = (
    0.069,
    0.0,
    0.069,
    0.0,
    0.01,
    0.0,
    0.0,
)

const alpha = (
    -pi/2,
    pi/2,
    -pi/2,
    pi/2,
    -pi/2,
    pi/2,
    0.0,
)

const Q = [
    0.0 -1.0 0.0 0.0
    1.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0
    0.0 0.0 0.0 0.0
]  # 偏微分演算行列

const g = [0.0, 0.0, -9.81, 0.0]'  # 重力加速度ベクトル（横ベクトル）

function r_bar(i::Int64)
    [
        x_bar[i]
        y_bar[i]
        z_bar[i]
        1.0
    ]
end

"""
(i-1)T(i)

thetaはq[i]です  
"""
function T(i::Int64, theta::Float64)
    [
        cos(theta) -cos(alpha[i])*sin(theta) sin(alpha[i])*sin(theta) a[i]*cos(theta)
        sin(theta) cos(alpha[i])*cos(theta) -sin(alpha[i])*cos(theta) a[i]*sin(theta)
        0.0 sin(alpha[i]) cos(alpha[i]) d[i]
        0.0 0.0 0.0 1.0
    ]
end


function T(i::Int64, j::Int64, q::Vector{TU}) where TU
    z = T(i, q[i])
    _z = similar(z)
    for k in i+1:j
        #z *= T(k, q[k])
        mul!(_z, T(k, q[k]), z)
        copy!(z, _z)
    end
    z
end


"""慣性モーメント"""
function J(i::Int64)
    [
        (-Ixx[i]+Iyy[i]+Izz[i])/2 Ixy[i]                   Ixz[i]                   m[i]*x_bar[i]
        Ixy[i]                    (Ixx[i]-Iyy[i]+Izz[i])/2 Iyz[i]                   m[i]*y_bar[i]
        Ixz[i]                    Iyz[i]                   (Ixx[i]+Iyy[i]-Izz[i])/2 m[i]*z_bar[i]
        m[i]*x_bar[i]             m[i]*y_bar[i]            m[i]*z_bar[i]            m[i]
    ]
end

# """慣性モーメント"""
# function J(i::Int64)
#     [
#         (-Ixx[i]+Iyy[i]+Izz[i])/2 Ixy[i] Ixz[i] m[i]*r_bars[1, i]
#         Ixy[i] (Ixx[i]-Iyy[i]+Izz[i])/2 Iyz[i] m[i]*r_bars[2, i]
#         Ixz[i] Iyz[i] (Ixx[i]+Iyy[i]-Izz[i])/2 m[i]*r_bars[3, i]
#         m[i]*r_bars[1,i] m[i]*r_bars[2,i] m[i]*r_bars[3,i] m[i]
#     ]
# end


function U(i::Int64, j::Int64, q::Vector{TU}) where TU
    if j <= i
        return T(1, j-1, q) * Q * T(j, i, q)
    else
        return zeros(TU, 4, 4)
    end
end


function U(i::Int64, j::Int64, k::Int64, q::Vector{TU}) where TU
    if i >= k >= j
        return T(1, j-1, q) * Q * T(j, k-1, q) * Q * T(k, i, q)
    elseif i >= j >= k
        return T(1, k-1, q) * Q * T(k, j-1, q) * Q * T(j, i, q)
    else
        return zeros(TU, 4, 4)
    end
end


### 慣性，コリオリ，重力 ###
"""慣性行列の要素"""
function M(i::Int64, k::Int64, q::Vector{TU}) where TU
    j_start = max(i, k)
    z = 0.0
    for j in j_start:n
        z += tr(U(j, k, q) * J(j) * U(j, i, q)')
    end
    z
end


function h(i::Int64, k::Int64, m::Int64, q::Vector{TU}) where TU
    j_start = max(i, k, m)
    z = 0.0
    for j in j_start:n
        z += tr(U(j, k, m, q) * J(i) * U(j, i, q)')
    end
    z
end


"""コリオリ行列の要素"""
function C(i::Int64, q::Vector{TU}, dq::Vector{TU}) where TU
    z = 0.0
    for k in 1:n
        for m = 1:n
            z += h(i, k, m, q) * dq[k] * dq[m]
        end
    end
    z
end


"""重力行列の要素"""
function G(i::Int64, q::Vector{TU}) where TU
    z = 0.0
    for j in i:n
        z += -m[j] * g * U(j, i, q) * r_bar(j)
    end
    z
end



### 目的の行列作成 ###

"""慣性行列"""
function M(q::Vector{TU}) where TU
    z = zeros(TU, 7, 7)
    for i in 1:7
        for j in 1:7
            z[i, j] = M(i, j, q)
        end
    end
    z
end

"""コリオリ行列"""
function C(q::Vector{TU}, dq::Vector{TU}) where TU
    z = zeros(TU, 7, 1)
    for i in 1:7
        z[i] = C(i, q, dq)
    end
    z
end

"""重力行列"""
function G(q::Vector{TU}) where TU
    z = zeros(TU, 7, 1)
    for i in 1:7
        z[i] = G(i, q)
    end
    z
end


"""トルクを計算

計算トルク法です  
q:（所望の）関節角度ベクトル  
dq:（所望の）関節角速度ベクトル  
ddq:（所望の）関節角速度ベクトル  
"""
function calc_torque(q::Vector{TU}, dq::Vector{TU}, ddq::Vector{TU}) where TU
    M(q)*ddq .+ C(q, dq) .+ G(q)
end


"""現実世界での加速度


u : 入力トルク  
F : 外力  
q:（現在の）関節角度ベクトル  
dq:（現在の）関節角速度ベクトル  
"""
function calc_real_ddq(u::Array{TU}, F::Array{TU}, q::Vector{TU}, dq::Vector{TU}) where TU
    inv(M(q)) * (u .+ F .- (C(q, dq) .+ G(q)))
end

#export calc_real_ddq

"""テスト用"""
function _test()
    q = [
        0.0
        -31
        0.0
        43.0
        0.0
        72.0
        0.0
    ] * pi / 180
    dq = q
    ddq = q

    u = calc_torque(q, dq, ddq)
    #println(u)
    F = zeros(Float64, 7, 1)
    r_ddq = calc_real_ddq(u, F, q, dq)
    #println(r_ddq)
    #println(typeof(r_ddq))
end

#end

@time for i in 1:100; _test(); end
#@time for i in 1:100; Dynamics._test(); end