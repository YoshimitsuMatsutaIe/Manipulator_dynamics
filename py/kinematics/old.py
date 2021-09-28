import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D 
from math import cos, sin, tan, pi
import math

class DHparam:
    """DHパラメータ"""
    
    def __init__(self, alpha, a, d, theta):
        self.alpha = alpha
        self.a = a
        self.d = d
        self.theta = theta


class HomogeneousTransformationMatrix:
    """同次変換行列"""
    
    def __init__(self, DHparam, M=None):
        self.update(DHparam, M)
        return
    
    def update(self, p, M):
        """情報を更新"""
        
        if M is not None:
            self.t = M  # 同次変換行列
        else:
            self.t = np.array([
                [cos(p.theta), -sin(p.theta), 0, p.a],
                [sin(p.theta)*cos(p.alpha), cos(p.theta)*cos(p.alpha), -sin(p.alpha), -p.d*sin(p.alpha)],
                [sin(p.theta)*sin(p.alpha), cos(p.theta)*sin(p.alpha), cos(p.alpha), p.d*cos(p.alpha)],
                [0, 0, 0, 1],
            ])  # 同次変換行列
        
        self.o = self.t[0:3, 3:4]
        self.R = self.t[0:3, 0:3]
        self.rx = self.t[0:3, 0:1]
        self.ry = self.t[0:3, 1:2]
        self.rz = self.t[0:3, 2:3]
        
        return
    
    def __mul__(self, other):
        M = self.t @ other.t
        return HomogeneousTransformationMatrix(DHparam=None, M=M)


# パラメータ
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

q_neutral = np.array([[0, -31, 0, 43, 0, 72, 0]]).T * pi/180  # ニュートラルの姿勢

qr = q_neutral  # 右手の関節角度ベクトル
ql = q_neutral  # 左手の関節角度ベクトル

DHparams_r = [
    DHparam(0, 0, 0, qr[0, 0]),
    DHparam(-pi/2, L1, 0, qr[1, 0]+pi/2),
    DHparam(pi/2, 0, L2, qr[2, 0]),
    DHparam(-pi/2, L3, 0, qr[3, 0]),
    DHparam(pi/2, 0, L4, qr[4, 0]),
    DHparam(-pi/2, L5, 0, qr[5, 0]),
    DHparam(pi/2, 0, 0, qr[6, 0]),
]

DHparams_l = [
    DHparam(0, 0, 0, ql[0, 0]),
    DHparam(-pi/2, L1, 0, ql[1, 0]+pi/2),
    DHparam(pi/2, 0, L2, ql[2, 0]),
    DHparam(-pi/2, L3, 0, ql[3, 0]),
    DHparam(pi/2, 0, L4, ql[4, 0]),
    DHparam(-pi/2, L5, 0, ql[5, 0]),
    DHparam(pi/2, 0, 0, ql[6, 0]),
]


# 制御点のローカル座標
r_bar_1 = [
    np.array([[0, L1/2, -L0/2, 1]]).T,
    np.array([[0, -L1/2, -L0/2, 1]]).T,
    np.array([[L1/2, 0, -L0/2, 1]]).T,
    np.array([[-L1/2, 0, -L0/2, 1]]).T,
]  # 1座標系からみた制御点位置

r_bar_2 = [
    np.array([[0, 0, L3/2, 1]]).T,
    np.array([[0, 0, -L3/2, 1]]).T,
]

r_bar_3 = [
    np.array([[0, L3/2, -L2*2/3, 1]]).T,
    np.array([[0, -L3/2, -L2*2/3, 1]]).T,
    np.array([[L3/2, 0, -L2*2/3, 1]]).T,
    np.array([[-L3/2, 0, -L2*2/3, 1]]).T,
    np.array([[0, L3/2, -L2*1/3, 1]]).T,
    np.array([[0, -L3/2, -L2*1/3, 1]]).T,
    np.array([[L3/2, 0, -L2*1/3, 1]]).T,
    np.array([[-L3/2, 0, -L2*1/3, 1]]).T,
]

r_bar_4 = [
    np.array([[0, 0, L3/2, 1]]).T,
    np.array([[0, 0, -L3/2, 1]]).T,
]

r_bar_5 = [
    np.array([[0, L5/2, -L4/2, 1]]).T,
    np.array([[0, -L5/2, -L4/2, 1]]).T,
    np.array([[L5/2, 0, -L4/2, 1]]).T,
    np.array([[-L5/2, 0, -L4/2, 1]]).T,
]

r_bar_6 = [
    np.array([[0, 0, L5/2, 1]]).T,
    np.array([[0, 0, -L5/2, 1]]).T,
]

r_bar_7 = [
    np.array([[0, L5/2, L6/2, 1]]).T,
    np.array([[0, -L5/2, L6/2, 1]]).T,
    np.array([[L5/2, 0, L6/2, 1]]).T,
    np.array([[-L5/2, 0, L6/2, 1]]).T,
]

r_bar_GL = [
    np.array([[0, 0, 0, 1]]).T
]

# 追加
r_bars = [
    r_bar_1, r_bar_2, r_bar_3, r_bar_4, r_bar_5, r_bar_6, r_bar_7, r_bar_GL,
]


# 同次変換行列
# 右手
T_BR_Wo = HomogeneousTransformationMatrix(
    DHparam=None,
    M=np.array([
        [-math.sqrt(2)/2, math.sqrt(2)/2, 0, -L,],
        [-math.sqrt(2)/2, -math.sqrt(2)/2, 0, -h,],
        [0, 0, 1, H,],
        [0, 0, 0, 1,],
    ])
)

T_0_BR = HomogeneousTransformationMatrix(
    DHparam=None,
    M=np.array([
        [1, 0, 0, 0,],
        [0, 1, 0, 0,],
        [0, 0, 1, L0,],
        [0, 0, 0, 1,],
    ])
)

# 左手
T_BL_Wo = HomogeneousTransformationMatrix(
    DHparam=None,
    M=np.array([
        [math.sqrt(2)/2, math.sqrt(2)/2, 0, L,],
        [-math.sqrt(2)/2, math.sqrt(2)/2, 0, -h,],
        [0, 0, 1, H,],
        [0, 0, 0, 1,],
    ])
)

T_0_BL = HomogeneousTransformationMatrix(
    DHparam=None,
    M=np.array([
        [1, 0, 0, 0,],
        [0, 1, 0, 0,],
        [0, 0, 1, L0,],
        [0, 0, 0, 1,],
    ])
)

T_GR_7 = HomogeneousTransformationMatrix(
    DHparam=None,
    M=np.array([
        [1, 0, 0, 0,],
        [0, 1, 0, 0,],
        [0, 0, 1, L6,],
        [0, 0, 0, 1,],
    ])
)


Trs = [T_BR_Wo, T_0_BR]
Tls = [T_BL_Wo, T_0_BL]


for param in DHparams_r:
    Trs.append(HomogeneousTransformationMatrix(DHparam=param))
for param in DHparams_l:
    Tls.append(HomogeneousTransformationMatrix(DHparam=param))

Trs.append(T_GR_7)
Tls.append(T_GR_7)


# Wo基準の同次変換行列を作成
Trs_ = []
Tls_ = []

for i, T in enumerate(Trs):
    if i == 0:
        Trs_.append(T)
    else:
        Trs_.append(Trs_[i-1] * T)
for i, T in enumerate(Tls):
    if i == 0:
        Tls_.append(T)
    else:
        Tls_.append(Tls_[i-1] * T)


# 試しに図示
# データ作成
xrs, yrs, zrs = [0], [0], [0]
for T in Trs:
    xrs.append(T.o[0, 0])
    yrs.append(T.o[1, 0])
    zrs.append(T.o[2, 0])

xls, yls, zls = [0], [0], [0]
for T in Tls:
    xls.append(T.o[0, 0])
    yls.append(T.o[1, 0])
    zls.append(T.o[2, 0])


fig = plt.figure()
ax = fig.gca(projection = '3d')
ax.grid(True)
ax.set_xlabel('X[m]')
ax.set_ylabel('Y[m]')
ax.set_zlabel('Z[m]')
ax.plot(xrs, yrs, zrs, ".-", label = "R-joints",)
ax.plot(xls, yls, zls, ".-", label = "L-joints",)
ax.legend()
ax.set_box_aspect((1,1,1))

plt.show()
