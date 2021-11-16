
import numpy as np
from math import pi, sin, cos, tan


class BaxterDynamics:
    """バクスターロボットの動力学"""
    
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
        0.0,
        0.3644,
        0.0,
        0.3743,
        0.0,
        0.2295,
    )

    a = (
        0.069,
        0.0,
        0.069,
        0.0,
        0.01,
        0.0,
        0.0,
    )

    alpha = (
        pi/2,
        pi/2,
        -pi/2,
        pi/2,
        -pi/2,
        pi/2,
        0.0,
    )

    Q = np.array([
        [0, -1, 0, 0,],
        [1, 0, 0, 0,],
        [0, 0, 0, 0,],
        [0, 0, 0, 0,],
    ])  # 偏微分演算行列

    g = np.array([[0, 0, -9.81, 0]])  # 重力加速度ベクトル（横ベクトル）

    def r_bar(self, i):
        i -= 1
        return np.array([
            [self.x_bar[i]],
            [self.y_bar[i]],
            [self.z_bar[i]],
            [1],
        ])


    def Ti(self, i, q):
        """
        同時変換行列 (i-1)T(i)
        
        thetaはq[i]です
        """
        i -= 1
        theta = q[i, 0]
        return np.array([
            [cos(theta), -cos(self.alpha[i])*sin(theta), sin(self.alpha[i])*sin(theta), self.a[i]*cos(theta)],
            [sin(theta), cos(self.alpha[i])*cos(theta), -sin(self.alpha[i])*cos(theta), self.a[i]*sin(theta)],
            [0, sin(self.alpha[i]), cos(self.alpha[i]), self.d[i]],
            [0, 0, 0, 1],
        ])


    def Tij(self, i, j, q):
        """同時変換行列"""
        z = self.Ti(i, q)
        for k in range(i+1, j+1):
            z = z @ self.Ti(k, q)
        return z


    def J(self, i):
        """慣性モーメント"""
        i -= 1
        return np.array([
            [(-self.Ixx[i]+self.Iyy[i]+self.Izz[i])/2, self.Ixy[i], self.Ixz[i], self.m[i]*self.x_bar[i],],
            [self.Ixy[i], (self.Ixx[i]-self.Iyy[i]+self.Izz[i])/2, self.Iyz[i], self.m[i]*self.y_bar[i],],
            [self.Ixz[i], self.Iyz[i], (self.Ixx[i]+self.Iyy[i]-self.Izz[i])/2, self.m[i]*self.z_bar[i],],
            [self.m[i]*self.x_bar[i], self.m[i]*self.y_bar[i], self.m[i]*self.z_bar[i], self.m[i],],
        ])


    def Uij(self, i, j, q):
        if j <= i:
            #print(self.Tij(1, j-1, q))
            return self.Tij(1, j-1, q) @ self.Q @ self.Tij(j, i, q)
        else:
            return np.zeros((4, 4))


    def Uijk(self, i, j, k, q):
        if i >= k and k >= j:
            return self.Tij(1, j-1, q) @ self.Q @ self.Tij(j, k-1, q) @ self.Q @ self.Tij(k, i, q)
        elif i >= j and j >= k:
            return self.Tij(1, k-1, q) @ self.Q @ self.Tij(k, j-1, q) @ self.Q @ self.Tij(j, i, q)
        else:
            return np.zeros((4, 4))



    ### 慣性，コリオリ，重力 ###
    def Mik(self, i, k, q):
        """慣性行列の要素"""
        j_start = max(i, k)
        z = 0
        for j in range(j_start, self.n+1):
            #print(j)
            #print(self.Uij(j, k, q))
            z += np.trace(self.Uij(j, k, q) @ self.J(j) @ self.Uij(j, i, q).T)
        
        return z


    def h(self, i, k, m, q):
        """?"""
        j_start = max(i, k, m)
        z = 0
        for j in range(j_start, self.n+1):
            z += np.trace(self.Uijk(j, k, m, q) @ self.J(i-1) @ self.Uij(j, i, q).T)
        
        return z


    def Ci(self, i, q, dq):
        """コリオリ行列の要素"""
        z = 0
        for k in range(1, self.n+1):
            for m in range(1, self.n+1):
                z += self.h(i, k, m, q) * dq[k-1, 0] * dq[m-1, 0]
        return z


    def Gi(self, i, q):
        """重力行列の要素"""
        z = 0
        for j in range(i, self.n+1):
            z += -self.m[j-1] * self.g @ self.Uij(j, i, q) @ self.r_bar(j)
        return z


    ### 目的の行列作成 ###
    def M(self, q):
        """慣性行列"""
        z = np.zeros((7, 7))
        for i in range(1, 7+1):
            for j in range(1, 7+1):
                z[i-1, j-1] = self.Mik(i, j, q)
        return z
    
    
    def C(self, q, dq):
        """コリオリ行列"""
        z = np.zeros((7, 1))
        for i in range(1, 7+1):
            z[i-1] = self.Ci(i, q, dq)
        
        return z
    
    
    def G(self, q):
        """重力行列"""
        z = np.zeros((7, 1))
        for i in range(1, 7+1):
            z[i-1] = self.Gi(i, q)
        
        return z


    def calc_torque(self, q, dq, ddq):
        """トルクを計算
        
        計算トルク法です  
        q : 関節角度ベクトル  
        dq : 関節角速度ベクトル  
        ddq : （所望の）関節角速度ベクトル  
        """
        return self.M(q) @ ddq + self.C(q, dq) + self.G(q)



    def calc_real_ddq(self, u, F, q, dq):
        """現実世界での加速度

        u : トルクベクトル R(7)  
        F : 外力ベクトル R(7)  
        q : （現在の）関節角度ベクトル  
        dq : （現在の）関節角速度ベクトル  
        """
        return np.linalg.inv(self.M(q)) @ (u + F - (self.C(q, dq) + self.G(q)))




"""テスト用"""
def _test():
    q = np.array([[0, -31, 0, 43, 0, 72, 0,]]).T * pi / 180
    dq = q
    ddq = q


    d = BaxterDynamics()
    
    #i = 2
    #j = 6
    #print(d.Tij(i, j, q))
    
    i = 3
    j = 5
    print(d.M(q))
    
    #print(d.M(q))
    
    u = d.calc_torque(q, dq, ddq)
    #print(u)
    #F = np.zeros((7, 1))
    #r_ddq = d.calc_real_ddq(u, F, q, dq)
    #println(r_ddq)
    #println(typeof(r_ddq))


if __name__ == "__main__":

    _test()

    #@time for i in 1:100; _test(); 
    #@time for i in 1:10; Dynamics._test(); 


    # using Profile
    # @profile for i in 1:10; Dynamics._test();