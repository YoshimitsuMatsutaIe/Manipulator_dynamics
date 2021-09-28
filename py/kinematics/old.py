import numpy as np
import matplotlib.pyplot as plt
from math import cos, sin, tan, pi

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
    DHparam(0, 0, 0, qr[0]),
    DHparam(-pi/2, L1, 0, qr[1]+pi/2),
    DHparam(pi/2, 0, L2, qr[2]),
    DHparam(-pi/2, L3, 0, qr[3]),
    DHparam(pi/2, 0, L4, qr[4]),
    DHparam(-pi/2, L5, 0, qr[5]),
    DHparam(pi/2, 0, 0, qr[6]),
]

DHparams_l = [
    DHparam(0, 0, 0, ql[0]),
    DHparam(-pi/2, L1, 0, ql[1]+pi/2),
    DHparam(pi/2, 0, L2, ql[2]),
    DHparam(-pi/2, L3, 0, ql[3]),
    DHparam(pi/2, 0, L4, ql[4]),
    DHparam(-pi/2, L5, 0, ql[5]),
    DHparam(pi/2, 0, 0, ql[6]),
]


if __name__ == "__main__":
    p = DHparam(1, 2, 3, 4)
    hoge = HomogeneousTransformationMatrix(p)
    print(hoge.T)
