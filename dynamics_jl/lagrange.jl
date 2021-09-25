"""ラグランジュ法による入力計算（テスト）
"""

using CPUTime
#using Symbolics
using LinearAlgebra

n = 7

Ixx = (
    0.0470910226,
    0.027885975,
    0.0266173355,
    0.0131822787,
    0.0166774282,
    0.0070053791,
    0.0008162135,
)

Iyy = (
    0.035959884,
    0.020787492,
    0.012480083,
    0.009268520,
    0.003746311,
    0.005527552,
    0.0008735012,
)

Izz = (
    0.0376697645,
    0.0117520941,
    0.0284435520,
    0.0071158268,
    0.0167545726,
    0.0038760715,
    0.0005494148,
)

Ixy = (
    -0.0061487003,
    -0.0001882199,
    -0.0039218988,
    -0.0001966341,
    -0.0001865762,
    0.0001534806,
    0.000128440,
)

Iyz = (
    -0.0007808689,
    0.0020767576,
    -0.001083893,
    0.000745949,
    0.0006473235,
    -0.0002111503,
    0.0001057726,
)

Ixz = (
    0.0001278755,
    -0.00030096397,
    0.0002927063,
    0.0003603617,
    0.0001840370,
    -0.0004438478,
    0.00018969891,
)

x_bar = (
    -0.05117,
    0.00269,
    -0.07176,
    0.00159,
    -0.01168,
    0.00697,
    0.005137,
)

y_bar = (
    0.07908,
    -0.00529,
    0.08149,
    -0.01117,
    0.13111,
    0.006,
    0.0009572,
)

z_bar = (
    0.00086,
    0.06845,
    0.00132,
    0.02618,
    0.0046,
    0.06048,
    -0.06682,
)

m = (
    5.70044,
    3.22698,
    4.31272,
    2.07206,
    2.24665,
    1.60979,
    0.54218,
)

d = (
    0.2703,
    0,
    0.3644,
    0,
    0.3743,
    0,
    0.2295,
)

a = (
    0.069,
    0,
    0.069,
    0,
    0.01,
    0,
    0,
)

alpha = (
    -pi/2,
    pi/2,
    -pi/2,
    pi/2,
    -pi/2,
    pi/2,
    0,
)

Q = [
    0 -1 0 0
    1 0 0 0
    0 0 0 0
    0 0 0 0
]  # 偏微分演算行列

g = [
    0
    0
    -9.81
    0
]'  # 重力加速度

#@variables q[1:7]
#@variables dq[1:7]

q = [
    0
    -31
    0
    43
    0
    72
    0
] * pi / 180

dq = [
    1
    2
    3
    4
    5
    6
    7
]

ddq = [
    1
    2
    3
    4
    5
    6
    7
]

function r_bar(i)
    [
        x_bar[i]
        y_bar[i]
        z_bar[i]
        1
    ]
end


function T(i)
    """(i-1)T(i)"""
    theta = q[i]
    [
        cos(theta) -cos(alpha[i])*sin(theta) sin(alpha[i])*sin(theta) a[i]*cos(theta)
        sin(theta) cos(alpha[i])*cos(theta) -sin(alpha[i])*cos(theta) a[i]*sin(theta)
        0 sin(alpha[i]) cos(alpha[i]) d[i]
        0 0 0 1
    ]
end

function T(i, j)
    z = T(i)
    for k in i+1:j
        z *= T(k)
    end
    z
end

#println(T(1, 5))


function J(i)
    [
        (-Ixx[i]+Iyy[i]+Izz[i])/2 Ixy[i] Ixz[i] m[i]*x_bar[i]
        Ixy[i] (Ixx[i]-Iyy[i]+Izz[i])/2 Iyz[i] m[i]*y_bar[i]
        Ixz[i] Iyz[i] (Ixx[i]+Iyy[i]-Izz[i])/2 m[i]*z_bar[i]
        m[i]*x_bar[i] m[i]*y_bar[i] m[i]*z_bar[i] m[i]
    ]
end


function U(i, j)
    if j <= i
        T(1, j-1) * Q * T(j, i)
    else
        zeros(4, 4)
    end
end

function U(i, j, k)
    if i >= k >= j
        T(1, j-1) * Q * T(j, k-1) * Q * T(k, i)
    elseif i >= j >= k
        T(1, k-1) * Q * T(k, j-1) * Q * T(j, i)
    else
        zeros(4, 4)
    end
end


### 慣性，コリオリ，重力 ###

function M(i, k)
    j_start = max(i, k)
    z = 0
    for j in j_start:n
        z += tr(U(j, k) * J(j) * U(j, i)')
    end
    z
end


function h(i, k, m)
    j_start = max(i, k, m)
    z = 0
    for j in j_start:n
        z += tr(U(j, k, m) * J(i) * U(j, i)')
    end
    z
end


function C(i)
    z = 0
    for k in 1:n
        for m = 1:n
            z += h(i, k, m) * dq[k] * dq[m]
        end
    end
    #simplify(z)
    z
end


function G(i)
    z = 0
    for j in i:n
        z += -m[j] * g * U(j, i) * r_bar(j)
    end
    z
end

#println(h(5, 7, 1))
#@time println(C(5))
#println(r_bar(4))

#hoge = C(5)



# ### モジュール化テスト ###
# open("test.txt", "w+") do f
#     println(f, string(C(5)))
# end


### 目的の行列作成 ###

# 慣性行列
function M()
    z = zeros(7, 7)
    for i in 1:7
        for j in 1:7
            z[i, j] = M(i, j)
        end
    end
    z
end

function C()
    z = zeros(7, 1)
    for i in 1:7
        z[i] = C(i)
    end
    z
end

function G()
    z = zeros(7, 1)
    for i in 1:7
        z[i] = G(i)
    end
    z
end


function calc_torque()
    M()*ddq + C() + G()
end

@time println(calc_torque())
#@time for i in 1:10; calc_torq(); end