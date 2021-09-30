"""RMPいろいろ"""


import numpy as np
import math
from math import pi, cos, sin, tan

#import rmp_fromGDS_attract_xi_M


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


class OriginalRMP:
    """論文[R1]のRMP
    論文そのまま
    """
    def __init__(self, **kwargs):
        # アトラクター加速度
        self.attract_max_speed = kwargs.pop('attract_max_speed')
        self.attract_gain = kwargs.pop('attract_gain')
        self.attrat_a_damp_r = kwargs.pop('attract_a_damp_r')
        # アトラクター計量
        self.attract_sigma_W = kwargs.pop('attract_sigma_W')
        self.attract_sigma_H = kwargs.pop('attract_sigma_H')
        self.attract_A_damp_r = kwargs.pop('attract_A_damp_r')
        
        # 障害物加速度
        self.obs_scale_rep = kwargs.pop('obs_scale_rep')
        self.obs_scale_damp = kwargs.pop('obs_scale_damp')
        self.obs_ratio = kwargs.pop('obs_ratio')
        self.obs_rep_gain = kwargs.pop('obs_rep_gain')
        # 障害物計量
        self.obs_r = kwargs.pop('obs_r')
        
        # ジョイント制限処理加速度
        self.jl_gamma_p = kwargs.pop('jl_gamma_p')
        self.jl_gamma_d = kwargs.pop('jl_gamma_d')
        # ジョイント制限処理計量
        self.jl_lambda = kwargs.pop('jl_lambda')
        self.joint_limit_upper = kwargs.pop('joint_limit_upper')
        self.joint_limit_lower = kwargs.pop('joint_limit_lower')
    
    
    ## アトラクターRMP
    def a_attract(self, z, dz, z0):
        """アトラクタ加速度"""
        
        max_speed = self.attract_max_speed
        gain = self.attract_gain
        damp = gain / max_speed
        alpha = self.attrat_a_damp_r
        
        return gain * soft_normal(z0 - z, alpha) - damp * dz
    
    def metric_attract(self, z, dz, z0, a):
        """アトラクタ計量"""
        
        sigma_W = self.attract_sigma_W
        sigma_H = self.attract_sigma_H
        alpha = self.attract_A_damp_r
        
        dis = np.linalg.norm(z0 - z)
        weight = np.exp(-dis / sigma_W)
        beta_attract = 1 - np.exp(-1 / 2 * (dis / sigma_H) ** 2)
        
        return weight * basic_metric_H(a, alpha, beta_attract) # 論文
    
    
    # 障害物RMP
    def a_obs(self, z, dz, z0):
        """障害物加速度"""
        
        scale_rep = self.obs_scale_rep  # 感知半径？
        scale_damp = self.obs_scale_damp
        ratio = self.obs_ratio
        rep_gain = self.obs_rep_gain
        damp_gain = rep_gain * ratio
        
        #x = z0 - z  # 反対かも？
        x = z - z0
        #dz = -dz  # 二乗するから関係なし
        
        dis = np.linalg.norm(x)  # 距離関数
        dis_grad = x / dis  # 距離関数の勾配
        
        # 斥力項．障害物に対する位置ベースの反発力？
        alpha_rep = rep_gain * np.exp(-dis / scale_rep)  # 斥力の活性化関数
        a_rep = alpha_rep * dis_grad  # 斥力項
        
        # ダンピング項．障害物に向かう速度にペナルティを課す？
        P_obs = max(0, -(dz).T @ dis_grad) * dis_grad @ dis_grad.T @ dz  # 零空間射影演算子
        alpha_damp = damp_gain / (dis / scale_damp + 1e-7)  # ダンピングの活性化関数
        a_damp = alpha_damp * P_obs  # ダンピング項
        #print("a_obs = ", a_rep + a_damp)
        return a_rep + a_damp
    
    def metric_obs(self, z, dz, z0, f_obs):
        """障害物計量"""
        
        r = self.obs_r  # この半径外に障害物が来ると発散
        x = z - z0
        dis = np.linalg.norm(x)
        weight_obs = (dis / r) ** 2 - 2 * dis / r + 1  #3次スプライン?
        #return weight_obs * basic_metric_H(f_obs, 1, 0.8)  # これかも？
        return weight_obs * np.eye(3, dtype = np.float32)  # 誤植？
        #return weight_obs * f_obs @ f_obs.T * 0.0001
    
    
    # ジョイント制限RMP
    def a_joint_limit(self, q, dq):
        """ジョイント制限処理加速度"""
        
        gamma_p = self.jl_gamma_p
        gamma_d = self.jl_gamma_d
        z = gamma_p * (-q) - gamma_d * dq
        #print("z = ", z)
        a = np.linalg.inv(
            D_sigma(
                q, 
                self.joint_limit_lower, 
                self.joint_limit_upper)) @ z
        return a
        
    
    def metric_joint_limit(self, q):
        """ジョイント制限処理計量"""
        dof = len(q)
        return self.jl_lambda * np.eye(dof)



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



# class RMPfromGDS:
#     """[R2],[R4]のRMP
    
#     注意
#     0 ＜ w_l ＜ w_u ＜ infty
#     alpha_fはたぶん減速半径
    
    
#     """
    
#     def __init__(self, **kwargs):
#         # アトラクターRMP関連
#         self.attract_max_speed = kwargs.pop('attract_max_speed')
#         self.attract_gain = kwargs.pop('attract_gain')
#         self.attract_alpha_f = kwargs.pop('attract_alpha_f')
#         self.attract_sigma_alpha = kwargs.pop('attract_sigma_alpha')
#         self.attract_sigma_gamma = kwargs.pop('attract_sigma_gamma')
#         self.attract_w_u = kwargs.pop('attract_w_u')
#         self.attract_w_l = kwargs.pop('attract_w_l')
#         self.attract_alpha = kwargs.pop('attract_alpha')
#         self.attract_epsilon = kwargs.pop('attract_epsilon')
        
#         # # 衝突回避関連
#         # self.obs_scale_rep = kwargs.pop('obs_scale_rep')
#         # self.obs_scale_damp = kwargs.pop('obs_scale_damp')
#         # self.obs_ratio = kwargs.pop('obs_ratio')
#         # self.obs_rep_gain = kwargs.pop('obs_rep_gain')
#         # # 障害物計量
#         # self.obs_r = kwargs.pop('obs_r')
        
#         # ジョイント制限処理力
#         self.jl_gamma_p = kwargs.pop('jl_gamma_p')
#         self.jl_gamma_d = kwargs.pop('jl_gamma_d')
#         # ジョイント制限処理計量
#         self.jl_lambda = kwargs.pop('jl_lambda')
#         self.joint_limit_upper = kwargs.pop('joint_limit_upper')
#         self.joint_limit_lower = kwargs.pop('joint_limit_lower')
#         self.jl_sigma = kwargs.pop('jl_sigma')
    
    
#     def inertia_attract(self, x, dx, x0, dx0):
#         """アトラクター慣性行列"""
#         z = x0 - x
#         dz = dx0 - dx
#         M = rmp_fromGDS_attract_xi_M.attract_M(
#             z, 
#             dz, 
#             self.attract_sigma_alpha, 
#             self.attract_sigma_gamma, 
#             self.attract_w_u, 
#             self.attract_w_l, 
#             self.attract_alpha, 
#             self.attract_epsilon)
#         return M
    
    
#     def f_attract(self, x, dx, x0, dx0, M_attract):
#         """アトラクト力（加速度？）"""
#         # パラメーター
#         gamma_p = self.attract_gain
#         gamma_d = self.attract_gain / self.attract_max_speed
#         alpha = self.attract_alpha_f
#         # 変数変換
#         z = x0 - x
#         dz = dx0 - dx
        
#         # メイン
#         f1 = -gamma_p * soft_normal(z, alpha) - gamma_d * dz
#         xi_M = rmp_fromGDS_attract_xi_M.attract_xi_M(
#             z, 
#             dz, 
#             self.attract_sigma_alpha, 
#             self.attract_sigma_gamma, 
#             self.attract_w_u, 
#             self.attract_w_l, 
#             self.attract_alpha, 
#             self.attract_epsilon)
#         carv = -np.linalg.inv(M_attract) @ xi_M
#         f = f1 + carv
#         return f
    
    
#     def metric_joint_limit(self, q, dq):
#         """ジョイント制限回避計量"""
#         A_ii = []
#         dof = len(q)
#         for i in range(0, dof, 1):
#             b = jl_b(
#                 q[i, 0], 
#                 dq[i, 0], 
#                 self.joint_limit_lower[i, 0], 
#                 self.joint_limit_upper[i, 0], 
#                 self.jl_sigma)
#             A_ii.append(b ** (-2))
#         A = self.jl_lambda * np.diag(A_ii)
#         return A
    
#     def f_joint_limit(self, q, dq, metric_jl):
#         """ジョイント制限処理力（加速度？）"""
        
#         gamma_p = self.jl_gamma_p
#         gamma_d = self.jl_gamma_d
#         ddq_old = gamma_p * (-q) - gamma_d * dq
        
#         A = metric_jl
#         dof = len(q)
        
#         d_xi_A = []
        
#         for i in range(0, dof, 1):
#             dAdq = dAiidqi(
#                 q[i, 0], 
#                 dq[i, 0], 
#                 self.joint_limit_lower[i, 0],
#                 self.joint_limit_upper[i, 0],
#                 self.jl_sigma)
#             d_xi_A.append(1/2 * dAdq * dq[i, 0])
        
#         xi_A = np.diag(d_xi_A)
        
#         carv = np.linalg.inv(A) @ xi_A
        
#         f = ddq_old - carv
#         #print("z = ", z)
#         return f




if __name__ == "__main__":
    print("Hello!!!")