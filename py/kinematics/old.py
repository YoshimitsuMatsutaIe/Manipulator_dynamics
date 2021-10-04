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
    
    def __init__(self, DHparam=None, M=None, zero=False):
        self.update(DHparam, M, zero)
        return
    
    def update(self, p, M, zero):
        """情報を更新"""
        
        if zero:
            self.t = np.zeros((4, 4))
        elif M is not None:
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


class BaxterKinematics:
    
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
    #q_neutral = np.array([[0, 31, 0, 43, 0, 72, 0]]).T * pi/180


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
    
    
    def __init__(self,):
        
        self.q_r = self.q_neutral  # 右手の関節角度ベクトル
        self.q_l = self.q_neutral  # 左手の関節角度ベクトル
        
        self.dq_r = self.q_neutral  # 右手の関節角速度ベクトル
        self.dq_l = self.q_neutral  # 左手の関節角速度ベクトル
        
        self.update_all()
        
    
    
    def update_all(self,):
        """全部アップデート"""
        
        self._update_HomogeneousTransformationMatrix(self.q_r, self.q_l)
        self._update_diff_HomogeneousTransformationMatrix()
        self._update_cpoints()
        self._update_jacobian()
        
        return
    
    
    def _update_HomogeneousTransformationMatrix(self, q_r, q_l):
        """同時変換行列を更新"""
    
        def update_DHparams(q_r, q_l):
            def DHparams(q):
                return [
                DHparam(0, 0, 0, q[0, 0]),
                DHparam(-pi/2, self.L1, 0, q[1, 0]+pi/2),
                DHparam(pi/2, 0, self.L2, q[2, 0]),
                DHparam(-pi/2, self.L3, 0, q[3, 0]),
                DHparam(pi/2, 0, self.L4, q[4, 0]),
                DHparam(-pi/2, self.L5, 0, q[5, 0]),
                DHparam(pi/2, 0, 0, q[6, 0]),
                ]
            
            self.DHparams_r = DHparams(q_r)
            self.DHparams_l = DHparams(q_l)
            return
        
    
        update_DHparams(q_r, q_l)
    
        # 同次変換行列（ローカル座標の）
        # 右手
        T_BR_Wo = HomogeneousTransformationMatrix(
            DHparam=None,
            M=np.array([
                [-math.sqrt(2)/2, math.sqrt(2)/2, 0, -self.L,],
                [-math.sqrt(2)/2, -math.sqrt(2)/2, 0, -self.h,],
                [0, 0, 1, self.H,],
                [0, 0, 0, 1,],
            ])
        )

        T_0_BR = HomogeneousTransformationMatrix(
            DHparam=None,
            M=np.array([
                [1, 0, 0, 0,],
                [0, 1, 0, 0,],
                [0, 0, 1, self.L0,],
                [0, 0, 0, 1,],
            ])
        )

        # 左手
        T_BL_Wo = HomogeneousTransformationMatrix(
            DHparam=None,
            M=np.array([
                [math.sqrt(2)/2, math.sqrt(2)/2, 0, self.L,],
                [-math.sqrt(2)/2, math.sqrt(2)/2, 0, -self.h,],
                [0, 0, 1, self.H,],
                [0, 0, 0, 1,],
            ])
        )

        T_0_BL = HomogeneousTransformationMatrix(
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


        self.Trs = [T_BR_Wo, T_0_BR]
        self.Tls = [T_BL_Wo, T_0_BL]


        for param in self.DHparams_r:
            self.Trs.append(HomogeneousTransformationMatrix(DHparam=param))
        for param in self.DHparams_l:
            self.Tls.append(HomogeneousTransformationMatrix(DHparam=param))

        self.Trs.append(T_GR_7)
        self.Tls.append(T_GR_7)


        # Wo基準の同次変換行列を作成
        self.Trs_Wo = []
        self.Tls_Wo = []

        for i, T in enumerate(self.Trs):
            if i == 0:
                self.Trs_Wo.append(T)
            else:
                self.Trs_Wo.append(self.Trs_Wo[-1] * T)
        for i, T in enumerate(self.Tls):
            if i == 0:
                self.Tls_Wo.append(T)
            else:
                self.Tls_Wo.append(self.Tls_Wo[-1] * T)


    def _update_diff_HomogeneousTransformationMatrix(self,):
        """微分同次変換行列？を更新 & ジョイントに関するヤコビ行列を作成"""
        
        def _update(Ts, Ts_Wo,):
            """4つのヤコビ行列を計算
            
            Ts : ローカル同時変換行列
            Ts_Wo : グローバル同次変換行列
            """
            
            dTj_dqis = []
            # j = 0
            for i in range(7):
                dTj_dqi = []
                for j in range(8):
                    
                    if j < i:
                        dTj_dqi.append(
                            HomogeneousTransformationMatrix(zero=True)
                        )
                    
                    elif j == i:
                        dTj_dqi.append(Ts_Wo[j+2] * self.A)
                    
                    else:
                        dTj_dqi.append(dTj_dqi[-1] * Ts[j+2])
                    
                dTj_dqis.append(dTj_dqi)
            
            
            Jaxs, Jays, Jazs, Jos = [], [], [], []
            for dT in [list(x) for x in zip(*dTj_dqis)]:
                _Jax = [T.rx_bar for T in dT]
                _Jay = [T.ry_bar for T in dT]
                _Jaz = [T.rz_bar for T in dT]
                _Jo = [T.o_bar for T in dT]
                
                Jax = np.concatenate(_Jax, axis=1)
                Jay = np.concatenate(_Jay, axis=1)
                Jaz = np.concatenate(_Jaz, axis=1)
                Jo = np.concatenate(_Jo, axis=1)

                Jaxs.append(Jax)
                Jays.append(Jay)
                Jazs.append(Jaz)
                Jos.append(Jo)
            
            return Jaxs, Jays, Jazs, Jos
        
        self.Jaxs_r, self.Jays_r, self.Jazs_r, self.Jos_r = _update(self.Trs, self.Trs_Wo)
        self.Jaxs_l, self.Jays_l, self.Jazs_l, self.Jos_l = _update(self.Tls, self.Tls_Wo)

        #print(self.Jaxs_l)
        
        return

    def _update_jacobian(self,):
        
        def _calc_Jo_global(Jax, Jay, Jaz, Jo, r_bar):
            z_bar = (Jax * r_bar[0,0] + Jay * r_bar[1,0] + Jaz * r_bar[2,0] + Jo)
            #print(z_bar)
            return z_bar[0:3, :]
        
        self.Jo_global_r = []
        for Jax, Jay, Jaz, Jo in zip(self.Jaxs_r, self.Jays_r, self.Jazs_r, self.Jos_r):
            self.Jo_global_r.append(_calc_Jo_global(Jax, Jay, Jaz, Jo, self.r_bar_zero))
        self.Jo_global_l = []
        for Jax, Jay, Jaz, Jo in zip(self.Jaxs_l, self.Jays_l, self.Jazs_l, self.Jos_l):
            self.Jo_global_l.append(_calc_Jo_global(Jax, Jay, Jaz, Jo, self.r_bar_zero))

        #print("Jo_global_l=", self.Jo_global_l)
        return


    def _update_cpoints(self,):
        
        # 制御点の位置を計算
        self.cpoints_r, self.cpoints_l = [], []
        for i, r_bar in enumerate(self.r_bars_all):
            c_temp = []
            T_temp = self.Trs_Wo[i+2]
            for r_bar_ in r_bar:
                c_temp.append(T_temp.t @ r_bar_)
            self.cpoints_r.append(c_temp)
        for i, r_bar in enumerate(self.r_bars_all):
            c_temp = []
            T_temp = self.Tls_Wo[i+2]
            for r_bar_ in r_bar:
                c_temp.append(T_temp.t @ r_bar_)
            self.cpoints_l.append(c_temp)
        return




    def get_joint_positions(self,):
        """ジョイント原点座標を取得"""
        return [
            [T.o for T in self.Trs_Wo],
            [T.o for T in self.Tls_Wo],
        ]

    def get_cpoint_positions(self,):
        """制御点座標を全取得"""
        return 





def main():
    start = time.time()
    kinema = BaxterKinematics()
    kinema.update_all()
    os = kinema.get_joint_positions()
    xrs, yrs, zrs = [0], [0], [0]
    for o in os[0]:
        xrs.append(o[0, 0])
        yrs.append(o[1, 0])
        zrs.append(o[2, 0])

    xls, yls, zls = [0], [0], [0]
    for o in os[1]:
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
    for i, cs in enumerate(kinema.cpoints_r):
        cs_ = np.concatenate(cs, axis=1)
        xs = cs_[0, :].tolist()
        ys = cs_[1, :].tolist()
        zs = cs_[2, :].tolist()
        ax.scatter(xs, ys, zs, label = "R-" + cs_name[i])
    for i, cs in enumerate(kinema.cpoints_l):
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