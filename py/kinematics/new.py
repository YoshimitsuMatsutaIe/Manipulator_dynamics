import numpy as np
import matplotlib.pyplot as plt
from math import cos, sin, tan, pi
import math
import time

from numpy.lib.shape_base import expand_dims

class DHparam:
    """DHパラメータ"""
    
    def __init__(self, alpha, a, d, theta):
        self.alpha = alpha
        self.a = a
        self.d = d
        self.theta = theta


class HomogeneousTransformationMatrix:
    """同次変換行列"""
    
    def __init__(self, DHparam=None, M=None,):
        self.update(DHparam, M,)
        return
    
    def update(self, p, M,):
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
        #self.R = self.t[0:3, 0:3]
        #self.rx = self.t[0:3, 0:1]
        #self.ry = self.t[0:3, 1:2]
        #self.rz = self.t[0:3, 2:3]
        self.rx_bar = self.t[:, 0:1]
        self.ry_bar = self.t[:, 1:2]
        self.rz_bar = self.t[:, 2:3]
        self.o_bar = self.t[:, 3:4]
        return

    def __mul__(self, other):
        M = self.t @ other.t
        return HomogeneousTransformationMatrix(DHparam=None, M=M)

    @classmethod
    def zero(cls,):
        return HomogeneousTransformationMatrix(
            M=np.zeros((4, 4))
        )


class BaxterRobotArmKinematics:
    
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


    # 制御点のローカル座標
    r_bars_in_1 = [
        np.array([[0, L1/2, -L0/2, 1]]).T,
        np.array([[0, -L1/2, -L0/2, 1]]).T,
        np.array([[L1/2, 0, -L0/2, 1]]).T,
        np.array([[-L1/2, 0, -L0/2, 1]]).T,
    ]  # 1座標系からみた制御点位置

    r_bars_in_2 = [
        np.array([[0, 0, L3/2, 1]]).T,
        np.array([[0, 0, -L3/2, 1]]).T,
    ]

    r_bars_in_3 = [
        np.array([[0, L3/2, -L2*2/3, 1]]).T,
        np.array([[0, -L3/2, -L2*2/3, 1]]).T,
        np.array([[L3/2, 0, -L2*2/3, 1]]).T,
        np.array([[-L3/2, 0, -L2*2/3, 1]]).T,
        np.array([[0, L3/2, -L2*1/3, 1]]).T,
        np.array([[0, -L3/2, -L2*1/3, 1]]).T,
        np.array([[L3/2, 0, -L2*1/3, 1]]).T,
        np.array([[-L3/2, 0, -L2*1/3, 1]]).T,
    ]

    r_bars_in_4 = [
        np.array([[0, 0, L3/2, 1]]).T,
        np.array([[0, 0, -L3/2, 1]]).T,
    ]

    r_bars_in_5 = [
        np.array([[0, L5/2, -L4/2, 1]]).T,
        np.array([[0, -L5/2, -L4/2, 1]]).T,
        np.array([[L5/2, 0, -L4/2, 1]]).T,
        np.array([[-L5/2, 0, -L4/2, 1]]).T,
    ]

    r_bars_in_6 = [
        np.array([[0, 0, L5/2, 1]]).T,
        np.array([[0, 0, -L5/2, 1]]).T,
    ]

    r_bars_in_7 = [
        np.array([[0, L5/2, L6/2, 1]]).T,
        np.array([[0, -L5/2, L6/2, 1]]).T,
        np.array([[L5/2, 0, L6/2, 1]]).T,
        np.array([[-L5/2, 0, L6/2, 1]]).T,
    ]

    r_bars_in_GL = [
        np.array([[0, 0, 0, 1]]).T
    ]

    # 追加
    r_bars_all = [
        r_bars_in_1, r_bars_in_2, r_bars_in_3, r_bars_in_4, r_bars_in_5, r_bars_in_6, r_bars_in_7, r_bars_in_GL,
    ]

    r_bar_zero = np.array([[0, 0, 0, 1]]).T
    

    A = HomogeneousTransformationMatrix(
        M=np.array([
            [0, -1, 0, 0],
            [1, 0, 0, 0],
            [0, 0, 0, 0],
            [0, 0, 0, 0],
        ])
    )  # 偏微分演算行列
    
    
    def __init__(self, isLeft):
        """
        
        isLeft : 左手か否か
        """
        
        self.isLeft = isLeft
        
        self.q = self.q_neutral  # 左手の関節角度ベクトル
        self.dq = self.q_neutral  # 左手の関節角速度ベクトル
        
        self.update_all()
        
    
    
    def update_all(self,):
        """全部アップデート"""

        self._update_HomogeneousTransformationMatrix(self.q)
        self._update_diff_HomogeneousTransformationMatrix()
        self._update_cpoints()
        self._update_jacobian()

        return
    
    
    def _update_HomogeneousTransformationMatrix(self, q):
        """同時変換行列を更新"""

        DHparams = [
            DHparam(0, 0, 0, q[0, 0]),
            DHparam(-pi/2, self.L1, 0, q[1, 0]+pi/2),
            DHparam(pi/2, 0, self.L2, q[2, 0]),
            DHparam(-pi/2, self.L3, 0, q[3, 0]),
            DHparam(pi/2, 0, self.L4, q[4, 0]),
            DHparam(-pi/2, self.L5, 0, q[5, 0]),
            DHparam(pi/2, 0, 0, q[6, 0]),
        ]

    
        # 同次変換行列（ローカル座標の）
        
        if self.isLeft:  # 左手
            T_BLorR_Wo = HomogeneousTransformationMatrix(
                DHparam=None,
                M=np.array([
                    [math.sqrt(2)/2, math.sqrt(2)/2, 0, self.L,],
                    [-math.sqrt(2)/2, math.sqrt(2)/2, 0, -self.h,],
                    [0, 0, 1, self.H,],
                    [0, 0, 0, 1,],
                ])
            )
            T_0_BLorR = HomogeneousTransformationMatrix(
                DHparam=None,
                M=np.array([
                    [1, 0, 0, 0,],
                    [0, 1, 0, 0,],
                    [0, 0, 1, self.L0,],
                    [0, 0, 0, 1,],
                ])
            )
        
        else:  #右手
            T_BLorR_Wo = HomogeneousTransformationMatrix(
                DHparam=None,
                M=np.array([
                    [-math.sqrt(2)/2, math.sqrt(2)/2, 0, -self.L,],
                    [-math.sqrt(2)/2, -math.sqrt(2)/2, 0, -self.h,],
                    [0, 0, 1, self.H,],
                    [0, 0, 0, 1,],
                ])
            )

            T_0_BLorR = HomogeneousTransformationMatrix(
                DHparam=None,
                M=np.array([
                    [1, 0, 0, 0,],
                    [0, 1, 0, 0,],
                    [0, 0, 1, self.L0,],
                    [0, 0, 0, 1,],
                ])
            )
        
        
        T_GR_7 = HomogeneousTransformationMatrix(
            DHparam=None,
            M=np.array([
                [1, 0, 0, 0,],
                [0, 1, 0, 0,],
                [0, 0, 1, self.L6,],
                [0, 0, 0, 1,],
            ])
        )


        self.Ts = [T_BLorR_Wo, T_0_BLorR]


        for param in DHparams:
            self.Ts.append(HomogeneousTransformationMatrix(DHparam=param))


        self.Ts.append(T_GR_7)


        # Wo基準の同次変換行列を作成
        self.Ts_Wo = []
        for i, T in enumerate(self.Ts):
            if i == 0:
                self.Ts_Wo.append(T)
            else:
                self.Ts_Wo.append(self.Ts_Wo[-1] * T)

        return


    def _update_diff_HomogeneousTransformationMatrix(self,):
        """微分同次変換行列？を更新 & ジョイントに関するヤコビ行列を作成"""

        dTj_dqis = []
        for i in range(7):
            dTj_dqi = []
            for j in range(8):
                
                if j < i:
                    dTj_dqi.append(
                        HomogeneousTransformationMatrix.zero()
                    )
                
                elif j == i:
                    dTj_dqi.append(self.Ts_Wo[j+2] * self.A)
                
                else:
                    dTj_dqi.append(dTj_dqi[-1] * self.Ts[j+2])
                
            dTj_dqis.append(dTj_dqi)
        
        
        self.Jaxs, self.Jays, self.Jazs, self.Jos = [], [], [], []
        for dT in [list(x) for x in zip(*dTj_dqis)]:
            _Jax = [T.rx_bar for T in dT]
            _Jay = [T.ry_bar for T in dT]
            _Jaz = [T.rz_bar for T in dT]
            _Jo = [T.o_bar for T in dT]
            
            Jax = np.concatenate(_Jax, axis=1)
            Jay = np.concatenate(_Jay, axis=1)
            Jaz = np.concatenate(_Jaz, axis=1)
            Jo = np.concatenate(_Jo, axis=1)

            self.Jaxs.append(Jax)
            self.Jays.append(Jay)
            self.Jazs.append(Jaz)
            self.Jos.append(Jo)
        
        return

    def _update_jacobian(self,):
        
        def _calc_Jo_global(Jax, Jay, Jaz, Jo, r_bar):
            z_bar = (Jax * r_bar[0,0] + Jay * r_bar[1,0] + Jaz * r_bar[2,0] + Jo)
            return z_bar[0:3, :]
        
        self.Jo_global = []
        for Jax, Jay, Jaz, Jo in zip(self.Jaxs, self.Jays, self.Jazs, self.Jos):
            self.Jo_global.append(_calc_Jo_global(Jax, Jay, Jaz, Jo, self.r_bar_zero))

        return


    def _update_cpoints(self,):
        
        # 制御点の位置を計算
        self.cpoints = []
        for i, r_bar in enumerate(self.r_bars_all):
            c_temp = []
            T_temp = self.Ts_Wo[i+2]
            for r_bar_ in r_bar:
                c_temp.append(T_temp.t @ r_bar_)
            self.cpoints.append(c_temp)
        return



    def get_joint_positions(self,):
        """ジョイント原点座標を取得"""
        return [T.o for T in self.Ts_Wo]

    def get_cpoint_positions(self,):
        """制御点座標を全取得"""
        return 





def main():
    start = time.time()
    right = BaxterRobotArmKinematics(isLeft=False)
    right.update_all()
    os_r = right.get_joint_positions()
    xrs, yrs, zrs = [0], [0], [0]
    for o in os_r:
        xrs.append(o[0, 0])
        yrs.append(o[1, 0])
        zrs.append(o[2, 0])

    left = BaxterRobotArmKinematics(isLeft=True)
    os_l = left.get_joint_positions()
    xls, yls, zls = [0], [0], [0]
    for o in os_l:
        xls.append(o[0, 0])
        yls.append(o[1, 0])
        zls.append(o[2, 0])



    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.grid(True)
    ax.set_xlabel('X[m]')
    ax.set_ylabel('Y[m]')
    ax.set_zlabel('Z[m]')
    ax.plot(xrs, yrs, zrs, ".-", label = "R-joints",)
    ax.plot(xls, yls, zls, ".-", label = "L-joints",)



    cs_name = ("1", "2", "3", "4", "5", "6", "7", "GL")
    for i, cs in enumerate(right.cpoints):
        cs_ = np.concatenate(cs, axis=1)
        xs = cs_[0, :].tolist()
        ys = cs_[1, :].tolist()
        zs = cs_[2, :].tolist()
        ax.scatter(xs, ys, zs, label = "R-" + cs_name[i])
    for i, cs in enumerate(left.cpoints):
        cs_ = np.concatenate(cs, axis=1)
        xs = cs_[0, :].tolist()
        ys = cs_[1, :].tolist()
        zs = cs_[2, :].tolist()
        ax.scatter(xs, ys, zs, label = "L-" + cs_name[i])


    ax.legend()
    ax.set_box_aspect((1, 1, 1))


    print("実行時間 ", time.time() - start)
    
    plt.show()


    return

if __name__ == "__main__":
    main()