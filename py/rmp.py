"""RMPいろいろ"""


import numpy as np
import math
from math import pi, cos, sin, tan

import rmp_fromGDS_attract_xi_M


class RMP:
    def __init__(self, goal_attractor, collision_avoidance,):
        self.goal_attractor = goal_attractor
        self.collision_avoidance = collision_avoidance



def pullback(f, M, J, dJ=None, dx= None):
    """pullback演算"""
    
    if dJ is None and dx is None:
        _f = J.T @ f
        _M = J.T @ M @ J
    else:

        _f = J.T @ (f - M @ dJ @ dx)
        _M = J.T @ M @ J
    return _f, _M



## マニピュレータの論文[R1]のやつ
def soft_normal(v, alpha):
    """ソフト正規化関数"""
    v_norm = np.linalg.norm(v)
    softmax = v_norm + 1 / alpha * np.log(1 + np.exp(-2 * alpha * v_norm))
    return v / softmax

def metric_stretch(v, alpha):
    """空間を一方向に伸ばす計量"""
    xi = soft_normal(v, alpha)
    return xi @ xi.T

def basic_metric_H(f, alpha, beta):
    """基本の計量"""
    f_norm = np.linalg.norm(f)
    f_softmax = f_norm + 1 / alpha * np.log(1 + np.exp(-2 * alpha * f_norm))
    s = f / f_softmax
    return beta * s @ s.T + (1 - beta) * np.eye(3)

def sigma_L(q, q_min, q_max):
    """アフィン変換andシグモイド変換写像"""
    return (q_max - q_min) * (1 / (1 + np.exp(-q))) + q_min

def D_sigma(q, q_min, q_max):
    """ジョイント制限に関する対角ヤコビ行列"""
    diags = (q_max - q_min) * (np.exp(-q) / (1 + np.exp(-q)) ** 2)
    #print(np.diag(diags.ravel()))
    return np.diag(diags.ravel())


# R1のやつ
class OriginalRMPAttractor:
    """論文[R1]のRMP"""

    def __init__(self, **kwargs):
        # アトラクター加速度
        self.attract_max_speed = kwargs.pop('attract_max_speed')
        self.attract_gain = kwargs.pop('attract_gain')
        self.attrat_a_damp_r = kwargs.pop('attract_a_damp_r')
        # アトラクター計量
        self.attract_sigma_W = kwargs.pop('attract_sigma_W')
        self.attract_sigma_H = kwargs.pop('attract_sigma_H')
        self.attract_A_damp_r = kwargs.pop('attract_A_damp_r')

    def _a(self, z, dz, z0):
        """アトラクタ加速度"""
        
        damp = self.attract_gain / self.attract_max_speed
        a = self.attract_gain * soft_normal(z0 - z, self.attrat_a_damp_r) - damp * dz
        return a
    
    def _metric(self, z, dz, z0, a):
        """アトラクタ計量"""
        
        dis = np.linalg.norm(z0 - z)
        weight = np.exp(-dis / self.attract_sigma_W)
        beta_attract = 1 - np.exp(-1 / 2 * (dis / self.attract_sigma_H) ** 2)
        
        return weight * basic_metric_H(a, self.attract_A_damp_r, beta_attract) # 論文

    def get_canonical(self, z, dz, z0):
        """form []"""
        
        a = self._a(z, dz, z0)
        M = self._metric(z, dz, z0, a)
        return a, M

    def get_natural(self, z, dz, z0, dz0):
        """form ()"""
        
        a, M = self.get_canonical(z, dz, z0)
        f = M @ a
        return f, M



class OriginalRMPCollisionAvoidance:
    """論文[R1]のRMP"""

    def __init__(self, **kwargs):
        # 障害物加速度
        self.obs_scale_rep = kwargs.pop('obs_scale_rep')  # 感知半径？
        self.obs_scale_damp = kwargs.pop('obs_scale_damp')
        self.obs_ratio = kwargs.pop('obs_ratio')
        self.obs_rep_gain = kwargs.pop('obs_rep_gain')
        # 障害物計量
        self.obs_r = kwargs.pop('obs_r')

    def _a(self, z, dz, z0):
        """障害物加速度"""
        
        damp_gain = self.obs_rep_gain * self.obs_ratio
        
        #x = z0 - z  # 反対かも？
        x = z - z0
        #dz = -dz  # 二乗するから関係なし
        
        dis = np.linalg.norm(x)  # 距離関数
        dis_grad = x / dis  # 距離関数の勾配
        
        # 斥力項．障害物に対する位置ベースの反発力？
        alpha_rep = self.obs_rep_gain * np.exp(-dis / self.obs_scale_rep)  # 斥力の活性化関数
        a_rep = alpha_rep * dis_grad  # 斥力項
        
        # ダンピング項．障害物に向かう速度にペナルティを課す？
        P_obs = max(0, -(dz).T @ dis_grad) * dis_grad @ dis_grad.T @ dz  # 零空間射影演算子
        alpha_damp = damp_gain / (dis / self.obs_scale_damp + 1e-7)  # ダンピングの活性化関数
        a_damp = alpha_damp * P_obs  # ダンピング項
        #print("a_obs = ", a_rep + a_damp)
        return a_rep + a_damp
    
    def _metric(self, z, dz, z0, f_obs):
        """障害物計量"""
        
        r = self.obs_r  # この半径外に障害物が来ると発散
        x = z - z0
        dis = np.linalg.norm(x)
        weight_obs = (dis / r) ** 2 - 2 * dis / r + 1  #3次スプライン?
        #return weight_obs * basic_metric_H(f_obs, 1, 0.8)  # これかも？
        return weight_obs * np.eye(3, dtype = np.float32)  # 誤植？
        #return weight_obs * f_obs @ f_obs.T * 0.0001

    def get_canonical(self, z, dz, z0, dz0):
        """form []"""
        
        a = self._a(z, dz, z0)
        M = self._metric(z, dz, z0, a)
        return a, M

    def get_natural(self, z, dz, z0, dz0):
        """form ()"""
        
        a, M = self.get_canonical(z, dz, z0, dz0)
        f = M @ a
        return f, M



class OriginalRMPJointLimitAvoidance:
    """論文[R1]のRMP"""

    def __init__(self, **kwargs):
        # ジョイント制限処理加速度
        self.jl_gamma_p = kwargs.pop('jl_gamma_p')
        self.jl_gamma_d = kwargs.pop('jl_gamma_d')
        # ジョイント制限処理計量
        self.jl_lambda = kwargs.pop('jl_lambda')

    def _a(self, q, dq, q_max, q_min):
        """ジョイント制限処理加速度"""
        
        z = self.jl_gamma_p * (-q) - self.jl_gamma_d * dq
        #print("z = ", z)
        a = np.linalg.inv(
            D_sigma(q, q_min, q_max)) @ z
        return a
    
    def _metric(self, q, dq, q_max, q_min,):
        """ジョイント制限処理計量"""
        dof = len(q)
        return self.jl_lambda * np.eye(dof)

    def get_canonical(self, q, dq, q_max, q_min):
        """form []"""
        
        a = self._a(q, dq, q_max, q_min)
        M = self._metric(q, dq, q_max, q_min,)
        return a, M

    def get_natural(self, q, dq, q_max, q_min):
        """form ()"""
        
        a, M = self.get_canonical(q, dq, q_max, q_min)
        f = M @ a
        return f, M




def jl_alpha_upper(dq, sigma):
    return 1 - math.exp(-(max(dq, 0) ** 2) / (2 * sigma ** 2))

def jl_alpha_lower(dq, sigma):
    return 1 - math.exp(-(min(dq, 0) ** 2) / (2 * sigma ** 2))

def jl_s(q, q_min, q_max):
    return (q - q_min) / (q_max - q_min)

def jl_b(q, dq, q_min, q_max, sigma):
    """???を計算"""
    s = jl_s(q, q_min, q_max)
    d = 4 * s * (1 - s)
    alpha_u = jl_alpha_upper(dq, sigma)
    alpha_l = jl_alpha_lower(dq, sigma)
    return s * (alpha_u * d + (1 - alpha_u)) + (1 - s) * (alpha_l * d + (1 - alpha_l))

def dAiidqi(q, dq, q_min, q_max, sigma):
    """????に使用"""
    alpha_l = jl_alpha_lower(dq, sigma)
    alpha_u = jl_alpha_upper(dq, sigma)
    z = 2*(q_max - q_min)**6*(-4*alpha_l*(q - q_max)*(q - q_min) - 4*alpha_l*(q - q_max)*(2*q - q_max - q_min) + 4*alpha_u*(q - q_max)*(q - q_min) + 4*alpha_u*(q - q_min)*(2*q - q_max - q_min)\
        - (1 - alpha_u)*(q_max - q_min)**2 - (alpha_l - 1)*(q_max - q_min)**2)\
            /((q - q_max)*(4*alpha_l*(q - q_max)*(q - q_min) + (alpha_l - 1)*(q_max - q_min)**2) - (q - q_min)*(4*alpha_u*(q - q_max)*(q - q_min) + (alpha_u - 1)*(q_max - q_min)**2))**3
    return z




class RMPfromGDSAttractor:
    """論文[R1]のRMP"""

    def __init__(self, **kwargs):
        self.attract_max_speed = kwargs.pop('attract_max_speed')
        self.attract_gain = kwargs.pop('attract_gain')
        self.attract_alpha_f = kwargs.pop('attract_alpha_f')
        self.attract_sigma_alpha = kwargs.pop('attract_sigma_alpha')
        self.attract_sigma_gamma = kwargs.pop('attract_sigma_gamma')
        self.attract_w_u = kwargs.pop('attract_w_u')
        self.attract_w_l = kwargs.pop('attract_w_l')
        self.attract_alpha = kwargs.pop('attract_alpha')
        self.attract_epsilon = kwargs.pop('attract_epsilon')

    def _inertia_matrix(self, x, dx, x0, dx0):
        """アトラクター慣性行列"""
        z = x0 - x
        dz = dx0 - dx
        M = rmp_fromGDS_attract_xi_M.attract_M(
            z, 
            dz, 
            self.attract_sigma_alpha, 
            self.attract_sigma_gamma, 
            self.attract_w_u, 
            self.attract_w_l, 
            self.attract_alpha, 
            self.attract_epsilon)
        return M
    
    def _f(self, x, dx, x0, dx0, M_attract):
        """アトラクト力（加速度？）"""
        # パラメーター
        gamma_p = self.attract_gain
        gamma_d = self.attract_gain / self.attract_max_speed
        alpha = self.attract_alpha_f
        # 変数変換
        z = x0 - x
        dz = dx0 - dx
        
        # メイン
        f1 = -gamma_p * soft_normal(z, alpha) - gamma_d * dz
        xi_M = rmp_fromGDS_attract_xi_M.attract_xi_M(
            z, 
            dz, 
            self.attract_sigma_alpha, 
            self.attract_sigma_gamma, 
            self.attract_w_u, 
            self.attract_w_l, 
            self.attract_alpha, 
            self.attract_epsilon
        )
        carv = -np.linalg.inv(M_attract) @ xi_M
        f = f1 + carv
        return f

    def get_natural(self, x, dx, x0, dx0):
        """form ()"""
        
        M = self._inertia_matrix(x, dx, x0, dx0)
        f = self._f(x, dx, x0, dx0, M)
        
        return f, M




class RMPfromGDSCollisionAvoidance:
    """論文[R2]のRMP"""

    def __init__(self, **kwargs):
        self.rw = kwargs.pop('rw')
        self.sigma = kwargs.pop('sigma')
        self.alpha = kwargs.pop('alpha')

    def _w(self, s):
        """重み関数"""
        # if s < self.rw:
        #     return (self.rw - s)**2 / s
        # else:
        #     return 0
        
        return 1/s**4

    def _dwds(self, s):
        """重み関数のs微分"""
        # if s < self.rw:
        #     return 1 - self.rw**2 / s**3
        # else:
        #     return 0

        return -4*s**-5
    
    def _u(self, ds):
        """速度依存計量の速度依存部分"""
        if ds < 0:
            return 1 - np.exp(-ds**2 / (2 * self.sigma**2))
            #return 0
        else:
            return 0
    
    def _dudsdot(self, ds,):
        if ds < 0:
            return -np.exp(-ds**2 / (2 * self.sigma**2)) * (-ds / self.sigma**2)
            #return 0
        else:
            return 0
    
    def _delta(self, s, ds,):
        return self._u(ds) + 1/2 * ds * self._dudsdot(ds)

    def _xi(self, s, ds):
        """曲率"""
        return 1/2 * self._u(ds) * self._dwds(s) * ds**2

    def _phi_1(self, s):
        """バリア型ポテンシャル（R2の2次元の例より）"""
        return -1/2 * self.alpha * s**-2

    def _grad_dphi_1(self, s):
        return 1 * self.alpha * s**-3

    def _inertia(self, s, ds):
        """障害物計量"""
        return self._w(s) * self._delta(s, ds)
    
    def _f(self, s, ds):
        """障害物力"""
        w = self._w(s)
        grad_phi = self._grad_dphi_1(s)
        xi = self._xi(s, ds)
        
        return -w * grad_phi - xi

    def get_natural(self, x, dx, x0, dx0):
        """form ()
        
        ・RMP-treeを実装後には直す
        """
        
        S = (x0 - x) 
        V = (dx0 - dx)
        s = np.linalg.norm(S)
        ds = 1/s * (S.T @ V)
        
        m = self._inertia(s, ds)
        f = self._f(s, ds,)
        
        J = -(x - dx).T / s
        dJ = -1 / s**2 * ((dx0 - dx).T - (x0 - dx).T*ds)
        
        f = J.T * (f - m * dJ @ dx)
        M = J.T @ m @ J
        
        # f = J.T * f
        # M = J.T @ m @ J
        
        # print('f = ', f)
        # print('M = ', M)
        return f, M




class RMPfromGDSJointLimitAvoidance:
    """論文[R1]のRMP"""

    def __init__(self, **kwargs):
        # ジョイント制限処理力
        self.jl_gamma_p = kwargs.pop('jl_gamma_p')
        self.jl_gamma_d = kwargs.pop('jl_gamma_d')
        # ジョイント制限処理計量
        self.jl_lambda = kwargs.pop('jl_lambda')
        self.jl_sigma = kwargs.pop('jl_sigma')

    def _inertia_matrix(self, q, dq, q_max, q_min):
        """ジョイント制限回避計量"""
        A_ii = []
        dof = len(q)
        for i in range(0, dof, 1):
            b = jl_b(
                q[i, 0], 
                dq[i, 0], 
                q_min[i, 0], 
                q_max[i, 0], 
                self.jl_sigma)
            A_ii.append(b ** (-2))
        A = self.jl_lambda * np.diag(A_ii)
        return A
    
    def _f(self, q, dq, q_max, q_min, metric_jl):
        """ジョイント制限処理力（加速度？）"""
        
        gamma_p = self.jl_gamma_p
        gamma_d = self.jl_gamma_d
        ddq_old = gamma_p * (-q) - gamma_d * dq
        
        A = metric_jl
        dof = len(q)
        
        d_xi_A = []
        
        for i in range(0, dof, 1):
            dAdq = dAiidqi(
                q[i, 0], 
                dq[i, 0], 
                q_min[i, 0],
                q_max[i, 0],
                self.jl_sigma)
            d_xi_A.append(1/2 * dAdq * dq[i, 0])
        
        xi_A = np.diag(d_xi_A)
        
        carv = np.linalg.inv(A) @ xi_A
        
        f = ddq_old - carv
        #print("z = ", z)
        return f

    def get_natural(self, q, dq, q_max, q_min):
        """form ()"""
        
        M = self._inertia_matrix(q, dq, q_max, q_min)
        f = self._f(q, dq, q_max, q_min, M)
        
        return f, M




if __name__ == "__main__":
    print("Hello!!!")