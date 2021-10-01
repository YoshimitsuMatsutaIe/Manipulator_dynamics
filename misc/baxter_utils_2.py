"""baxterロボットの色々関数とクラス
・順運動学写像，ヤコビ行列など
"""

import numpy as np
from math import cos, sin, tan, pi, sqrt


class BaxterKinematics2:
    """ローカル座標系の原点，ヤコビ行列等を計算"""
    
    def __init__(self, **kwargs):
        self.L = kwargs.pop('L')
        self.h = kwargs.pop('h')
        self.H = kwargs.pop('H')
        self.L0 = kwargs.pop('L0')
        self.L1 = kwargs.pop('L1')
        self.L2 = kwargs.pop('L2')
        self.L3 = kwargs.pop('L3')
        self.L4 = kwargs.pop('L4')
        self.L5 = kwargs.pop('L5')
        self.L6 = kwargs.pop('L6')
    
    def origin_Wo(self, q):
        L, h, H, L0, L1, L2, L3, L4, L5, L6 = self.L, self.h, self.H, self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6
        q1, q2, q3, q4, q5, q6, q7 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0], q[6, 0]
        z = np.zeros((3, 1))
        return z
    
    def origin_BL(self, q):
        L, h, H, L0, L1, L2, L3, L4, L5, L6 = self.L, self.h, self.H, self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6
        q1, q2, q3, q4, q5, q6, q7 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0], q[6, 0]
        z = np.array([[L], 
                      [-h], 
                      [H]])
        return z
    
    def origin_0(self, q):
        L, h, H, L0, L1, L2, L3, L4, L5, L6 = self.L, self.h, self.H, self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6
        q1, q2, q3, q4, q5, q6, q7 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0], q[6, 0]
        z = np.array([[L], 
                      [-h], 
                      [H + L0]])
        return z
    
    def origin_1(self, q):
        L, h, H, L0, L1, L2, L3, L4, L5, L6 = self.L, self.h, self.H, self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6
        q1, q2, q3, q4, q5, q6, q7 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0], q[6, 0]
        z = np.array([[L], 
                      [-h], 
                      [H + L0]])
        return z
    
    def origin_2(self, q):
        L, h, H, L0, L1, L2, L3, L4, L5, L6 = self.L, self.h, self.H, self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6
        q1, q2, q3, q4, q5, q6, q7 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0], q[6, 0]
        z = np.array([[L + L1*cos(q1 + pi/4)], 
                      [L1*sin(q1 + pi/4) - h], 
                      [H + L0]])
        return z
    
    def origin_3(self, q):
        L, h, H, L0, L1, L2, L3, L4, L5, L6 = self.L, self.h, self.H, self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6
        q1, q2, q3, q4, q5, q6, q7 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0], q[6, 0]
        z = np.array([[L + L1*cos(q1 + pi/4) + L2*cos(q2)*cos(q1 + pi/4)], 
                      [L1*sin(q1 + pi/4) + L2*sin(q1 + pi/4)*cos(q2) - h], 
                      [H + L0 - L2*sin(q2)]])
        return z
    
    def origin_4(self, q):
        L, h, H, L0, L1, L2, L3, L4, L5, L6 = self.L, self.h, self.H, self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6
        q1, q2, q3, q4, q5, q6, q7 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0], q[6, 0]
        z = np.array([[L + L1*cos(q1 + pi/4) + L2*cos(q2)*cos(q1 + pi/4) - L3*sin(q2)*cos(q3)*cos(q1 + pi/4) - L3*sin(q3)*sin(q1 + pi/4)], 
                      [L1*sin(q1 + pi/4) + L2*sin(q1 + pi/4)*cos(q2) - L3*sin(q2)*sin(q1 + pi/4)*cos(q3) + L3*sin(q3)*cos(q1 + pi/4) - h], 
                      [H + L0 - L2*sin(q2) - L3*cos(q2)*cos(q3)]])  
        return z
    
    def origin_5(self, q):
        L, h, H, L0, L1, L2, L3, L4, L5, L6 = self.L, self.h, self.H, self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6
        q1, q2, q3, q4, q5, q6, q7 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0], q[6, 0]
        z = np.array([[L + L1*cos(q1 + pi/4) + L2*cos(q2)*cos(q1 + pi/4) - L3*sin(q2)*cos(q3)*cos(q1 + pi/4) - L3*sin(q3)*sin(q1 + pi/4) - L4*sin(q2)*sin(q4)*cos(q3)*cos(q1 + pi/4) - L4*sin(q3)*sin(q4)*sin(q1 + pi/4) + L4*cos(q2)*cos(q4)*cos(q1 + pi/4)], 
                      [L1*sin(q1 + pi/4) + L2*sin(q1 + pi/4)*cos(q2) - L3*sin(q2)*sin(q1 + pi/4)*cos(q3) + L3*sin(q3)*cos(q1 + pi/4) - L4*sin(q2)*sin(q4)*sin(q1 + pi/4)*cos(q3) + L4*sin(q3)*sin(q4)*cos(q1 + pi/4) + L4*sin(q1 + pi/4)*cos(q2)*cos(q4) - h], 
                      [H + L0 - L2*sin(q2) - L3*cos(q2)*cos(q3) - L4*(sin(q2)*cos(q4) + sin(q4)*cos(q2)*cos(q3))]])
        return z
    
    def origin_6(self, q):
        L, h, H, L0, L1, L2, L3, L4, L5, L6 = self.L, self.h, self.H, self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6
        q1, q2, q3, q4, q5, q6, q7 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0], q[6, 0]
        z = np.array([[L + L1*cos(q1 + pi/4) + L2*cos(q2)*cos(q1 + pi/4) - L3*sin(q2)*cos(q3)*cos(q1 + pi/4) - L3*sin(q3)*sin(q1 + pi/4) - L4*sin(q2)*sin(q4)*cos(q3)*cos(q1 + pi/4) - L4*sin(q3)*sin(q4)*sin(q1 + pi/4) + L4*cos(q2)*cos(q4)*cos(q1 + pi/4) + L5*sin(q2)*sin(q3)*sin(q5)*cos(q1 + pi/4) - L5*sin(q2)*cos(q3)*cos(q4)*cos(q5)*cos(q1 + pi/4) - L5*sin(q3)*sin(q1 + pi/4)*cos(q4)*cos(q5) - L5*sin(q4)*cos(q2)*cos(q5)*cos(q1 + pi/4) - L5*sin(q5)*sin(q1 + pi/4)*cos(q3)], 
                      [L1*sin(q1 + pi/4) + L2*sin(q1 + pi/4)*cos(q2) - L3*sin(q2)*sin(q1 + pi/4)*cos(q3) + L3*sin(q3)*cos(q1 + pi/4) - L4*sin(q2)*sin(q4)*sin(q1 + pi/4)*cos(q3) + L4*sin(q3)*sin(q4)*cos(q1 + pi/4) + L4*sin(q1 + pi/4)*cos(q2)*cos(q4) + L5*sin(q2)*sin(q3)*sin(q5)*sin(q1 + pi/4) - L5*sin(q2)*sin(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) + L5*sin(q3)*cos(q4)*cos(q5)*cos(q1 + pi/4) - L5*sin(q4)*sin(q1 + pi/4)*cos(q2)*cos(q5) + L5*sin(q5)*cos(q3)*cos(q1 + pi/4) - h], 
                      [H + L0 - L2*sin(q2) - L3*cos(q2)*cos(q3) - L4*(sin(q2)*cos(q4) + sin(q4)*cos(q2)*cos(q3)) + L5*((sin(q2)*sin(q4) - cos(q2)*cos(q3)*cos(q4))*cos(q5) + sin(q3)*sin(q5)*cos(q2))]])
        return z
    
    def origin_7(self, q):
        L, h, H, L0, L1, L2, L3, L4, L5, L6 = self.L, self.h, self.H, self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6
        q1, q2, q3, q4, q5, q6, q7 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0], q[6, 0]
        z = np.array([[L + L1*cos(q1 + pi/4) + L2*cos(q2)*cos(q1 + pi/4) - L3*sin(q2)*cos(q3)*cos(q1 + pi/4) - L3*sin(q3)*sin(q1 + pi/4) - L4*sin(q2)*sin(q4)*cos(q3)*cos(q1 + pi/4) - L4*sin(q3)*sin(q4)*sin(q1 + pi/4) + L4*cos(q2)*cos(q4)*cos(q1 + pi/4) + L5*sin(q2)*sin(q3)*sin(q5)*cos(q1 + pi/4) - L5*sin(q2)*cos(q3)*cos(q4)*cos(q5)*cos(q1 + pi/4) - L5*sin(q3)*sin(q1 + pi/4)*cos(q4)*cos(q5) - L5*sin(q4)*cos(q2)*cos(q5)*cos(q1 + pi/4) - L5*sin(q5)*sin(q1 + pi/4)*cos(q3)], 
                      [L1*sin(q1 + pi/4) + L2*sin(q1 + pi/4)*cos(q2) - L3*sin(q2)*sin(q1 + pi/4)*cos(q3) + L3*sin(q3)*cos(q1 + pi/4) - L4*sin(q2)*sin(q4)*sin(q1 + pi/4)*cos(q3) + L4*sin(q3)*sin(q4)*cos(q1 + pi/4) + L4*sin(q1 + pi/4)*cos(q2)*cos(q4) + L5*sin(q2)*sin(q3)*sin(q5)*sin(q1 + pi/4) - L5*sin(q2)*sin(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) + L5*sin(q3)*cos(q4)*cos(q5)*cos(q1 + pi/4) - L5*sin(q4)*sin(q1 + pi/4)*cos(q2)*cos(q5) + L5*sin(q5)*cos(q3)*cos(q1 + pi/4) - h], 
                      [H + L0 - L2*sin(q2) - L3*cos(q2)*cos(q3) - L4*(sin(q2)*cos(q4) + sin(q4)*cos(q2)*cos(q3)) + L5*((sin(q2)*sin(q4) - cos(q2)*cos(q3)*cos(q4))*cos(q5) + sin(q3)*sin(q5)*cos(q2))]])
        return z
    
    def origin_GL(self, q):
        L, h, H, L0, L1, L2, L3, L4, L5, L6 = self.L, self.h, self.H, self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6
        q1, q2, q3, q4, q5, q6, q7 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0], q[6, 0]
        z = np.array([[L + L1*cos(q1 + pi/4) + L2*cos(q2)*cos(q1 + pi/4) - L3*sin(q2)*cos(q3)*cos(q1 + pi/4) - L3*sin(q3)*sin(q1 + pi/4) - L4*sin(q2)*sin(q4)*cos(q3)*cos(q1 + pi/4) - L4*sin(q3)*sin(q4)*sin(q1 + pi/4) + L4*cos(q2)*cos(q4)*cos(q1 + pi/4) + L5*sin(q2)*sin(q3)*sin(q5)*cos(q1 + pi/4) - L5*sin(q2)*cos(q3)*cos(q4)*cos(q5)*cos(q1 + pi/4) - L5*sin(q3)*sin(q1 + pi/4)*cos(q4)*cos(q5) - L5*sin(q4)*cos(q2)*cos(q5)*cos(q1 + pi/4) - L5*sin(q5)*sin(q1 + pi/4)*cos(q3) + L6*sin(q2)*sin(q3)*sin(q5)*sin(q6)*cos(q1 + pi/4) - L6*sin(q2)*sin(q4)*cos(q3)*cos(q6)*cos(q1 + pi/4) - L6*sin(q2)*sin(q6)*cos(q3)*cos(q4)*cos(q5)*cos(q1 + pi/4) - L6*sin(q3)*sin(q4)*sin(q1 + pi/4)*cos(q6) - L6*sin(q3)*sin(q6)*sin(q1 + pi/4)*cos(q4)*cos(q5) - L6*sin(q4)*sin(q6)*cos(q2)*cos(q5)*cos(q1 + pi/4) - L6*sin(q5)*sin(q6)*sin(q1 + pi/4)*cos(q3) + L6*cos(q2)*cos(q4)*cos(q6)*cos(q1 + pi/4)], 
                      [L1*sin(q1 + pi/4) + L2*sin(q1 + pi/4)*cos(q2) - L3*sin(q2)*sin(q1 + pi/4)*cos(q3) + L3*sin(q3)*cos(q1 + pi/4) - L4*sin(q2)*sin(q4)*sin(q1 + pi/4)*cos(q3) + L4*sin(q3)*sin(q4)*cos(q1 + pi/4) + L4*sin(q1 + pi/4)*cos(q2)*cos(q4) + L5*sin(q2)*sin(q3)*sin(q5)*sin(q1 + pi/4) - L5*sin(q2)*sin(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) + L5*sin(q3)*cos(q4)*cos(q5)*cos(q1 + pi/4) - L5*sin(q4)*sin(q1 + pi/4)*cos(q2)*cos(q5) + L5*sin(q5)*cos(q3)*cos(q1 + pi/4) + L6*sin(q2)*sin(q3)*sin(q5)*sin(q6)*sin(q1 + pi/4) - L6*sin(q2)*sin(q4)*sin(q1 + pi/4)*cos(q3)*cos(q6) - L6*sin(q2)*sin(q6)*sin(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) + L6*sin(q3)*sin(q4)*cos(q6)*cos(q1 + pi/4) + L6*sin(q3)*sin(q6)*cos(q4)*cos(q5)*cos(q1 + pi/4) - L6*sin(q4)*sin(q6)*sin(q1 + pi/4)*cos(q2)*cos(q5) + L6*sin(q5)*sin(q6)*cos(q3)*cos(q1 + pi/4) + L6*sin(q1 + pi/4)*cos(q2)*cos(q4)*cos(q6) - h], 
                      [H + L0 - L2*sin(q2) - L3*cos(q2)*cos(q3) - L4*(sin(q2)*cos(q4) + sin(q4)*cos(q2)*cos(q3)) + L5*((sin(q2)*sin(q4) - cos(q2)*cos(q3)*cos(q4))*cos(q5) + sin(q3)*sin(q5)*cos(q2)) + L6*(((sin(q2)*sin(q4) - cos(q2)*cos(q3)*cos(q4))*cos(q5) + sin(q3)*sin(q5)*cos(q2))*sin(q6) - (sin(q2)*cos(q4) + sin(q4)*cos(q2)*cos(q3))*cos(q6))]])
        return z
    
    def origins(self, q):
        o_Wo = self.origin_Wo(q)
        o_BL = self.origin_BL(q)
        o_0 = self.origin_0(q)
        o_1 = self.origin_1(q)
        o_2 = self.origin_2(q)
        o_3 = self.origin_3(q)
        o_4 = self.origin_4(q)
        o_5 = self.origin_5(q)
        o_6 = self.origin_6(q)
        o_7 = self.origin_7(q)
        o_GL = self.origin_GL(q)
        return [o_Wo, o_BL, o_0, o_1, o_2, o_3, o_4, o_5, o_6, o_7, o_GL]
    
    
    def jacobi_Wo(self, q):
        L, h, H, L0, L1, L2, L3, L4, L5, L6 = self.L, self.h, self.H, self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6
        q1, q2, q3, q4, q5, q6, q7 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0], q[6, 0]
        z = np.array([[0, 0, 0, 0, 0, 0, 0], 
                      [0, 0, 0, 0, 0, 0, 0], 
                      [0, 0, 0, 0, 0, 0, 0]])
        return z
    
    def jacobi_BL(self, q):
        L, h, H, L0, L1, L2, L3, L4, L5, L6 = self.L, self.h, self.H, self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6
        q1, q2, q3, q4, q5, q6, q7 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0], q[6, 0]
        z = np.array([[0, 0, 0, 0, 0, 0, 0], 
                      [0, 0, 0, 0, 0, 0, 0], 
                      [0, 0, 0, 0, 0, 0, 0]])
        return z
    
    def jacobi_0(self, q):
        L, h, H, L0, L1, L2, L3, L4, L5, L6 = self.L, self.h, self.H, self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6
        q1, q2, q3, q4, q5, q6, q7 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0], q[6, 0]
        z = np.array([[0, 0, 0, 0, 0, 0, 0], 
                      [0, 0, 0, 0, 0, 0, 0], 
                      [0, 0, 0, 0, 0, 0, 0]])
        return z
    
    def jacobi_1(self, q):
        L, h, H, L0, L1, L2, L3, L4, L5, L6 = self.L, self.h, self.H, self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6
        q1, q2, q3, q4, q5, q6, q7 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0], q[6, 0]
        z = np.array([[0, 0, 0, 0, 0, 0, 0], 
                      [0, 0, 0, 0, 0, 0, 0], 
                      [0, 0, 0, 0, 0, 0, 0]])
        return z
    
    def jacobi_2(self, q):
        L, h, H, L0, L1, L2, L3, L4, L5, L6 = self.L, self.h, self.H, self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6
        q1, q2, q3, q4, q5, q6, q7 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0], q[6, 0]
        z = np.array([[-L1*sin(q1 + pi/4), 0, 0, 0, 0, 0, 0], 
                      [L1*cos(q1 + pi/4), 0, 0, 0, 0, 0, 0], 
                      [0, 0, 0, 0, 0, 0, 0]])
        return z
    
    def jacobi_3(self, q):
        L, h, H, L0, L1, L2, L3, L4, L5, L6 = self.L, self.h, self.H, self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6
        q1, q2, q3, q4, q5, q6, q7 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0], q[6, 0]
        z = np.array([[-(L1 + L2*cos(q2))*sin(q1 + pi/4), -L2*sin(q2)*cos(q1 + pi/4), 0, 0, 0, 0, 0], 
                      [(L1 + L2*cos(q2))*cos(q1 + pi/4), -L2*sin(q2)*sin(q1 + pi/4), 0, 0, 0, 0, 0], 
                      [0, -L2*cos(q2), 0, 0, 0, 0, 0]])
        return z
    
    def jacobi_4(self, q):
        L, h, H, L0, L1, L2, L3, L4, L5, L6 = self.L, self.h, self.H, self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6
        q1, q2, q3, q4, q5, q6, q7 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0], q[6, 0]
        z = np.array([[-L1*sin(q1 + pi/4) - L2*sin(q1 + pi/4)*cos(q2) + L3*sin(q2)*sin(q1 + pi/4)*cos(q3) - L3*sin(q3)*cos(q1 + pi/4), -(L2*sin(q2) + L3*cos(q2)*cos(q3))*cos(q1 + pi/4), L3*(sin(q2)*sin(q3)*cos(q1 + pi/4) - sin(q1 + pi/4)*cos(q3)), 0, 0, 0, 0], 
                      [L1*cos(q1 + pi/4) + L2*cos(q2)*cos(q1 + pi/4) - L3*sin(q2)*cos(q3)*cos(q1 + pi/4) - L3*sin(q3)*sin(q1 + pi/4), -(L2*sin(q2) + L3*cos(q2)*cos(q3))*sin(q1 + pi/4), L3*(sin(q2)*sin(q3)*sin(q1 + pi/4) + cos(q3)*cos(q1 + pi/4)), 0, 0, 0, 0], 
                      [0, -L2*cos(q2) + L3*sin(q2)*cos(q3), L3*sin(q3)*cos(q2), 0, 0, 0, 0]])
        return z
    
    def jacobi_5(self, q):
        L, h, H, L0, L1, L2, L3, L4, L5, L6 = self.L, self.h, self.H, self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6
        q1, q2, q3, q4, q5, q6, q7 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0], q[6, 0]
        z = np.array([[-L1*sin(q1 + pi/4) - L2*sin(q1 + pi/4)*cos(q2) + L3*sin(q2)*sin(q1 + pi/4)*cos(q3) - L3*sin(q3)*cos(q1 + pi/4) + L4*sin(q2)*sin(q4)*sin(q1 + pi/4)*cos(q3) - L4*sin(q3)*sin(q4)*cos(q1 + pi/4) - L4*sin(q1 + pi/4)*cos(q2)*cos(q4), 
                       -(L2*sin(q2) + L3*cos(q2)*cos(q3) + L4*sin(q2)*cos(q4) + L4*sin(q4)*cos(q2)*cos(q3))*cos(q1 + pi/4), 
                       L3*sin(q2)*sin(q3)*cos(q1 + pi/4) - L3*sin(q1 + pi/4)*cos(q3) + L4*sin(q2)*sin(q3)*sin(q4)*cos(q1 + pi/4) - L4*sin(q4)*sin(q1 + pi/4)*cos(q3), 
                       -L4*(sin(q2)*cos(q3)*cos(q4)*cos(q1 + pi/4) + sin(q3)*sin(q1 + pi/4)*cos(q4) + sin(q4)*cos(q2)*cos(q1 + pi/4)), 0, 0, 0], 
                      [L1*cos(q1 + pi/4) + L2*cos(q2)*cos(q1 + pi/4) - L3*sin(q2)*cos(q3)*cos(q1 + pi/4) - L3*sin(q3)*sin(q1 + pi/4) - L4*sin(q2)*sin(q4)*cos(q3)*cos(q1 + pi/4) - L4*sin(q3)*sin(q4)*sin(q1 + pi/4) + L4*cos(q2)*cos(q4)*cos(q1 + pi/4), 
                       -(L2*sin(q2) + L3*cos(q2)*cos(q3) + L4*sin(q2)*cos(q4) + L4*sin(q4)*cos(q2)*cos(q3))*sin(q1 + pi/4), 
                       L3*sin(q2)*sin(q3)*sin(q1 + pi/4) + L3*cos(q3)*cos(q1 + pi/4) + L4*sin(q2)*sin(q3)*sin(q4)*sin(q1 + pi/4) + L4*sin(q4)*cos(q3)*cos(q1 + pi/4), 
                       L4*(-sin(q2)*sin(q1 + pi/4)*cos(q3)*cos(q4) + sin(q3)*cos(q4)*cos(q1 + pi/4) - sin(q4)*sin(q1 + pi/4)*cos(q2)), 0, 0, 0], 
                      [0, -L2*cos(q2) + L3*sin(q2)*cos(q3) + L4*(sin(q2)*sin(q4)*cos(q3) - cos(q2)*cos(q4)), (L3 + L4*sin(q4))*sin(q3)*cos(q2), L4*(sin(q2)*sin(q4) - cos(q2)*cos(q3)*cos(q4)), 0, 0, 0]])
        return z
    
    def jacobi_6(self, q):
        L, h, H, L0, L1, L2, L3, L4, L5, L6 = self.L, self.h, self.H, self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6
        q1, q2, q3, q4, q5, q6, q7 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0], q[6, 0]
        z = np.array([[-L1*sin(q1 + pi/4) - L2*sin(q1 + pi/4)*cos(q2) + L3*sin(q2)*sin(q1 + pi/4)*cos(q3) - L3*sin(q3)*cos(q1 + pi/4) + L4*sin(q2)*sin(q4)*sin(q1 + pi/4)*cos(q3) - L4*sin(q3)*sin(q4)*cos(q1 + pi/4) - L4*sin(q1 + pi/4)*cos(q2)*cos(q4) - L5*sin(q2)*sin(q3)*sin(q5)*sin(q1 + pi/4) + L5*sin(q2)*sin(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) - L5*sin(q3)*cos(q4)*cos(q5)*cos(q1 + pi/4) + L5*sin(q4)*sin(q1 + pi/4)*cos(q2)*cos(q5) - L5*sin(q5)*cos(q3)*cos(q1 + pi/4), 
                       (-L2*sin(q2) - L3*cos(q2)*cos(q3) - L4*sin(q2)*cos(q4) - L4*sin(q4)*cos(q2)*cos(q3) + L5*sin(q2)*sin(q4)*cos(q5) + L5*sin(q3)*sin(q5)*cos(q2) - L5*cos(q2)*cos(q3)*cos(q4)*cos(q5))*cos(q1 + pi/4), 
                       L3*sin(q2)*sin(q3)*cos(q1 + pi/4) - L3*sin(q1 + pi/4)*cos(q3) + L4*sin(q2)*sin(q3)*sin(q4)*cos(q1 + pi/4) - L4*sin(q4)*sin(q1 + pi/4)*cos(q3) + L5*sin(q2)*sin(q3)*cos(q4)*cos(q5)*cos(q1 + pi/4) + L5*sin(q2)*sin(q5)*cos(q3)*cos(q1 + pi/4) + L5*sin(q3)*sin(q5)*sin(q1 + pi/4) - L5*sin(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5), 
                       -L4*sin(q2)*cos(q3)*cos(q4)*cos(q1 + pi/4) - L4*sin(q3)*sin(q1 + pi/4)*cos(q4) - L4*sin(q4)*cos(q2)*cos(q1 + pi/4) + L5*sin(q2)*sin(q4)*cos(q3)*cos(q5)*cos(q1 + pi/4) + L5*sin(q3)*sin(q4)*sin(q1 + pi/4)*cos(q5) - L5*cos(q2)*cos(q4)*cos(q5)*cos(q1 + pi/4), 
                       L5*(sin(q2)*sin(q3)*cos(q5)*cos(q1 + pi/4) + sin(q2)*sin(q5)*cos(q3)*cos(q4)*cos(q1 + pi/4) + sin(q3)*sin(q5)*sin(q1 + pi/4)*cos(q4) + sin(q4)*sin(q5)*cos(q2)*cos(q1 + pi/4) - sin(q1 + pi/4)*cos(q3)*cos(q5)), 
                       0, 
                       0], 
                      [L1*cos(q1 + pi/4) + L2*cos(q2)*cos(q1 + pi/4) - L3*sin(q2)*cos(q3)*cos(q1 + pi/4) - L3*sin(q3)*sin(q1 + pi/4) - L4*sin(q2)*sin(q4)*cos(q3)*cos(q1 + pi/4) - L4*sin(q3)*sin(q4)*sin(q1 + pi/4) + L4*cos(q2)*cos(q4)*cos(q1 + pi/4) + L5*sin(q2)*sin(q3)*sin(q5)*cos(q1 + pi/4) - L5*sin(q2)*cos(q3)*cos(q4)*cos(q5)*cos(q1 + pi/4) - L5*sin(q3)*sin(q1 + pi/4)*cos(q4)*cos(q5) - L5*sin(q4)*cos(q2)*cos(q5)*cos(q1 + pi/4) - L5*sin(q5)*sin(q1 + pi/4)*cos(q3), 
                       (-L2*sin(q2) - L3*cos(q2)*cos(q3) - L4*sin(q2)*cos(q4) - L4*sin(q4)*cos(q2)*cos(q3) + L5*sin(q2)*sin(q4)*cos(q5) + L5*sin(q3)*sin(q5)*cos(q2) - L5*cos(q2)*cos(q3)*cos(q4)*cos(q5))*sin(q1 + pi/4), 
                       L3*sin(q2)*sin(q3)*sin(q1 + pi/4) + L3*cos(q3)*cos(q1 + pi/4) + L4*sin(q2)*sin(q3)*sin(q4)*sin(q1 + pi/4) + L4*sin(q4)*cos(q3)*cos(q1 + pi/4) + L5*sin(q2)*sin(q3)*sin(q1 + pi/4)*cos(q4)*cos(q5) + L5*sin(q2)*sin(q5)*sin(q1 + pi/4)*cos(q3) - L5*sin(q3)*sin(q5)*cos(q1 + pi/4) + L5*cos(q3)*cos(q4)*cos(q5)*cos(q1 + pi/4), 
                       -L4*sin(q2)*sin(q1 + pi/4)*cos(q3)*cos(q4) + L4*sin(q3)*cos(q4)*cos(q1 + pi/4) - L4*sin(q4)*sin(q1 + pi/4)*cos(q2) + L5*sin(q2)*sin(q4)*sin(q1 + pi/4)*cos(q3)*cos(q5) - L5*sin(q3)*sin(q4)*cos(q5)*cos(q1 + pi/4) - L5*sin(q1 + pi/4)*cos(q2)*cos(q4)*cos(q5), 
                       L5*(sin(q2)*sin(q3)*sin(q1 + pi/4)*cos(q5) + sin(q2)*sin(q5)*sin(q1 + pi/4)*cos(q3)*cos(q4) - sin(q3)*sin(q5)*cos(q4)*cos(q1 + pi/4) + sin(q4)*sin(q5)*sin(q1 + pi/4)*cos(q2) + cos(q3)*cos(q5)*cos(q1 + pi/4)), 
                       0, 
                       0], 
                      [0, 
                       -L2*cos(q2) + L3*sin(q2)*cos(q3) + L4*(sin(q2)*sin(q4)*cos(q3) - cos(q2)*cos(q4)) + L5*((sin(q2)*cos(q3)*cos(q4) + sin(q4)*cos(q2))*cos(q5) - sin(q2)*sin(q3)*sin(q5)), 
                       (L3*sin(q3) + L4*sin(q3)*sin(q4) + L5*(sin(q3)*cos(q4)*cos(q5) + sin(q5)*cos(q3)))*cos(q2), 
                       L4*(sin(q2)*sin(q4) - cos(q2)*cos(q3)*cos(q4)) + L5*(sin(q2)*cos(q4) + sin(q4)*cos(q2)*cos(q3))*cos(q5), 
                       -L5*((sin(q2)*sin(q4) - cos(q2)*cos(q3)*cos(q4))*sin(q5) - sin(q3)*cos(q2)*cos(q5)), 0, 0]])
        return z
    
    def jacobi_7(self, q):
        L, h, H, L0, L1, L2, L3, L4, L5, L6 = self.L, self.h, self.H, self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6
        q1, q2, q3, q4, q5, q6, q7 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0], q[6, 0]
        z = np.array([[-L1*sin(q1 + pi/4) - L2*sin(q1 + pi/4)*cos(q2) + L3*sin(q2)*sin(q1 + pi/4)*cos(q3) - L3*sin(q3)*cos(q1 + pi/4) + L4*sin(q2)*sin(q4)*sin(q1 + pi/4)*cos(q3) - L4*sin(q3)*sin(q4)*cos(q1 + pi/4) - L4*sin(q1 + pi/4)*cos(q2)*cos(q4) - L5*sin(q2)*sin(q3)*sin(q5)*sin(q1 + pi/4) + L5*sin(q2)*sin(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) - L5*sin(q3)*cos(q4)*cos(q5)*cos(q1 + pi/4) + L5*sin(q4)*sin(q1 + pi/4)*cos(q2)*cos(q5) - L5*sin(q5)*cos(q3)*cos(q1 + pi/4), 
                       (-L2*sin(q2) - L3*cos(q2)*cos(q3) - L4*sin(q2)*cos(q4) - L4*sin(q4)*cos(q2)*cos(q3) + L5*sin(q2)*sin(q4)*cos(q5) + L5*sin(q3)*sin(q5)*cos(q2) - L5*cos(q2)*cos(q3)*cos(q4)*cos(q5))*cos(q1 + pi/4), 
                       L3*sin(q2)*sin(q3)*cos(q1 + pi/4) - L3*sin(q1 + pi/4)*cos(q3) + L4*sin(q2)*sin(q3)*sin(q4)*cos(q1 + pi/4) - L4*sin(q4)*sin(q1 + pi/4)*cos(q3) + L5*sin(q2)*sin(q3)*cos(q4)*cos(q5)*cos(q1 + pi/4) + L5*sin(q2)*sin(q5)*cos(q3)*cos(q1 + pi/4) + L5*sin(q3)*sin(q5)*sin(q1 + pi/4) - L5*sin(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5), 
                       -L4*sin(q2)*cos(q3)*cos(q4)*cos(q1 + pi/4) - L4*sin(q3)*sin(q1 + pi/4)*cos(q4) - L4*sin(q4)*cos(q2)*cos(q1 + pi/4) + L5*sin(q2)*sin(q4)*cos(q3)*cos(q5)*cos(q1 + pi/4) + L5*sin(q3)*sin(q4)*sin(q1 + pi/4)*cos(q5) - L5*cos(q2)*cos(q4)*cos(q5)*cos(q1 + pi/4), 
                       L5*(sin(q2)*sin(q3)*cos(q5)*cos(q1 + pi/4) + sin(q2)*sin(q5)*cos(q3)*cos(q4)*cos(q1 + pi/4) + sin(q3)*sin(q5)*sin(q1 + pi/4)*cos(q4) + sin(q4)*sin(q5)*cos(q2)*cos(q1 + pi/4) - sin(q1 + pi/4)*cos(q3)*cos(q5)), 
                       0, 
                       0], 
                      [L1*cos(q1 + pi/4) + L2*cos(q2)*cos(q1 + pi/4) - L3*sin(q2)*cos(q3)*cos(q1 + pi/4) - L3*sin(q3)*sin(q1 + pi/4) - L4*sin(q2)*sin(q4)*cos(q3)*cos(q1 + pi/4) - L4*sin(q3)*sin(q4)*sin(q1 + pi/4) + L4*cos(q2)*cos(q4)*cos(q1 + pi/4) + L5*sin(q2)*sin(q3)*sin(q5)*cos(q1 + pi/4) - L5*sin(q2)*cos(q3)*cos(q4)*cos(q5)*cos(q1 + pi/4) - L5*sin(q3)*sin(q1 + pi/4)*cos(q4)*cos(q5) - L5*sin(q4)*cos(q2)*cos(q5)*cos(q1 + pi/4) - L5*sin(q5)*sin(q1 + pi/4)*cos(q3), 
                       (-L2*sin(q2) - L3*cos(q2)*cos(q3) - L4*sin(q2)*cos(q4) - L4*sin(q4)*cos(q2)*cos(q3) + L5*sin(q2)*sin(q4)*cos(q5) + L5*sin(q3)*sin(q5)*cos(q2) - L5*cos(q2)*cos(q3)*cos(q4)*cos(q5))*sin(q1 + pi/4), 
                       L3*sin(q2)*sin(q3)*sin(q1 + pi/4) + L3*cos(q3)*cos(q1 + pi/4) + L4*sin(q2)*sin(q3)*sin(q4)*sin(q1 + pi/4) + L4*sin(q4)*cos(q3)*cos(q1 + pi/4) + L5*sin(q2)*sin(q3)*sin(q1 + pi/4)*cos(q4)*cos(q5) + L5*sin(q2)*sin(q5)*sin(q1 + pi/4)*cos(q3) - L5*sin(q3)*sin(q5)*cos(q1 + pi/4) + L5*cos(q3)*cos(q4)*cos(q5)*cos(q1 + pi/4), 
                       -L4*sin(q2)*sin(q1 + pi/4)*cos(q3)*cos(q4) + L4*sin(q3)*cos(q4)*cos(q1 + pi/4) - L4*sin(q4)*sin(q1 + pi/4)*cos(q2) + L5*sin(q2)*sin(q4)*sin(q1 + pi/4)*cos(q3)*cos(q5) - L5*sin(q3)*sin(q4)*cos(q5)*cos(q1 + pi/4) - L5*sin(q1 + pi/4)*cos(q2)*cos(q4)*cos(q5), 
                       L5*(sin(q2)*sin(q3)*sin(q1 + pi/4)*cos(q5) + sin(q2)*sin(q5)*sin(q1 + pi/4)*cos(q3)*cos(q4) - sin(q3)*sin(q5)*cos(q4)*cos(q1 + pi/4) + sin(q4)*sin(q5)*sin(q1 + pi/4)*cos(q2) + cos(q3)*cos(q5)*cos(q1 + pi/4)), 
                       0, 
                       0], 
                      [0, 
                       -L2*cos(q2) + L3*sin(q2)*cos(q3) + L4*(sin(q2)*sin(q4)*cos(q3) - cos(q2)*cos(q4)) + L5*((sin(q2)*cos(q3)*cos(q4) + sin(q4)*cos(q2))*cos(q5) - sin(q2)*sin(q3)*sin(q5)), 
                       (L3*sin(q3) + L4*sin(q3)*sin(q4) + L5*(sin(q3)*cos(q4)*cos(q5) + sin(q5)*cos(q3)))*cos(q2), 
                       L4*(sin(q2)*sin(q4) - cos(q2)*cos(q3)*cos(q4)) + L5*(sin(q2)*cos(q4) + sin(q4)*cos(q2)*cos(q3))*cos(q5), 
                       -L5*((sin(q2)*sin(q4) - cos(q2)*cos(q3)*cos(q4))*sin(q5) - sin(q3)*cos(q2)*cos(q5)), 
                       0, 
                       0]])
        return z
    
    def jacobi_GL(self, q):
        L, h, H, L0, L1, L2, L3, L4, L5, L6 = self.L, self.h, self.H, self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6
        q1, q2, q3, q4, q5, q6, q7 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0], q[6, 0]
        z = np.array([[-L1*sin(q1 + pi/4) - L2*sin(q1 + pi/4)*cos(q2) + L3*sin(q2)*sin(q1 + pi/4)*cos(q3) - L3*sin(q3)*cos(q1 + pi/4) + L4*sin(q2)*sin(q4)*sin(q1 + pi/4)*cos(q3) - L4*sin(q3)*sin(q4)*cos(q1 + pi/4) - L4*sin(q1 + pi/4)*cos(q2)*cos(q4) - L5*sin(q2)*sin(q3)*sin(q5)*sin(q1 + pi/4) + L5*sin(q2)*sin(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) - L5*sin(q3)*cos(q4)*cos(q5)*cos(q1 + pi/4) + L5*sin(q4)*sin(q1 + pi/4)*cos(q2)*cos(q5) - L5*sin(q5)*cos(q3)*cos(q1 + pi/4) - L6*sin(q2)*sin(q3)*sin(q5)*sin(q6)*sin(q1 + pi/4) + L6*sin(q2)*sin(q4)*sin(q1 + pi/4)*cos(q3)*cos(q6) + L6*sin(q2)*sin(q6)*sin(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) - L6*sin(q3)*sin(q4)*cos(q6)*cos(q1 + pi/4) - L6*sin(q3)*sin(q6)*cos(q4)*cos(q5)*cos(q1 + pi/4) + L6*sin(q4)*sin(q6)*sin(q1 + pi/4)*cos(q2)*cos(q5) - L6*sin(q5)*sin(q6)*cos(q3)*cos(q1 + pi/4) - L6*sin(q1 + pi/4)*cos(q2)*cos(q4)*cos(q6), 
                       (-L2*sin(q2) - L3*cos(q2)*cos(q3) - L4*sin(q2)*cos(q4) - L4*sin(q4)*cos(q2)*cos(q3) + L5*sin(q2)*sin(q4)*cos(q5) + L5*sin(q3)*sin(q5)*cos(q2) - L5*cos(q2)*cos(q3)*cos(q4)*cos(q5) + L6*sin(q2)*sin(q4)*sin(q6)*cos(q5) - L6*sin(q2)*cos(q4)*cos(q6) + L6*sin(q3)*sin(q5)*sin(q6)*cos(q2) - L6*sin(q4)*cos(q2)*cos(q3)*cos(q6) - L6*sin(q6)*cos(q2)*cos(q3)*cos(q4)*cos(q5))*cos(q1 + pi/4), 
                       L3*sin(q2)*sin(q3)*cos(q1 + pi/4) - L3*sin(q1 + pi/4)*cos(q3) + L4*sin(q2)*sin(q3)*sin(q4)*cos(q1 + pi/4) - L4*sin(q4)*sin(q1 + pi/4)*cos(q3) + L5*sin(q2)*sin(q3)*cos(q4)*cos(q5)*cos(q1 + pi/4) + L5*sin(q2)*sin(q5)*cos(q3)*cos(q1 + pi/4) + L5*sin(q3)*sin(q5)*sin(q1 + pi/4) - L5*sin(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) + L6*sin(q2)*sin(q3)*sin(q4)*cos(q6)*cos(q1 + pi/4) + L6*sin(q2)*sin(q3)*sin(q6)*cos(q4)*cos(q5)*cos(q1 + pi/4) + L6*sin(q2)*sin(q5)*sin(q6)*cos(q3)*cos(q1 + pi/4) + L6*sin(q3)*sin(q5)*sin(q6)*sin(q1 + pi/4) - L6*sin(q4)*sin(q1 + pi/4)*cos(q3)*cos(q6) - L6*sin(q6)*sin(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5), 
                       -L4*sin(q2)*cos(q3)*cos(q4)*cos(q1 + pi/4) - L4*sin(q3)*sin(q1 + pi/4)*cos(q4) - L4*sin(q4)*cos(q2)*cos(q1 + pi/4) + L5*sin(q2)*sin(q4)*cos(q3)*cos(q5)*cos(q1 + pi/4) + L5*sin(q3)*sin(q4)*sin(q1 + pi/4)*cos(q5) - L5*cos(q2)*cos(q4)*cos(q5)*cos(q1 + pi/4) + L6*sin(q2)*sin(q4)*sin(q6)*cos(q3)*cos(q5)*cos(q1 + pi/4) - L6*sin(q2)*cos(q3)*cos(q4)*cos(q6)*cos(q1 + pi/4) + L6*sin(q3)*sin(q4)*sin(q6)*sin(q1 + pi/4)*cos(q5) - L6*sin(q3)*sin(q1 + pi/4)*cos(q4)*cos(q6) - L6*sin(q4)*cos(q2)*cos(q6)*cos(q1 + pi/4) - L6*sin(q6)*cos(q2)*cos(q4)*cos(q5)*cos(q1 + pi/4), 
                       L5*sin(q2)*sin(q3)*cos(q5)*cos(q1 + pi/4) + L5*sin(q2)*sin(q5)*cos(q3)*cos(q4)*cos(q1 + pi/4) + L5*sin(q3)*sin(q5)*sin(q1 + pi/4)*cos(q4) + L5*sin(q4)*sin(q5)*cos(q2)*cos(q1 + pi/4) - L5*sin(q1 + pi/4)*cos(q3)*cos(q5) + L6*sin(q2)*sin(q3)*sin(q6)*cos(q5)*cos(q1 + pi/4) + L6*sin(q2)*sin(q5)*sin(q6)*cos(q3)*cos(q4)*cos(q1 + pi/4) + L6*sin(q3)*sin(q5)*sin(q6)*sin(q1 + pi/4)*cos(q4) + L6*sin(q4)*sin(q5)*sin(q6)*cos(q2)*cos(q1 + pi/4) - L6*sin(q6)*sin(q1 + pi/4)*cos(q3)*cos(q5), 
                       L6*(sin(q2)*sin(q3)*sin(q5)*cos(q6)*cos(q1 + pi/4) + sin(q2)*sin(q4)*sin(q6)*cos(q3)*cos(q1 + pi/4) - sin(q2)*cos(q3)*cos(q4)*cos(q5)*cos(q6)*cos(q1 + pi/4) + sin(q3)*sin(q4)*sin(q6)*sin(q1 + pi/4) - sin(q3)*sin(q1 + pi/4)*cos(q4)*cos(q5)*cos(q6) - sin(q4)*cos(q2)*cos(q5)*cos(q6)*cos(q1 + pi/4) - sin(q5)*sin(q1 + pi/4)*cos(q3)*cos(q6) - sin(q6)*cos(q2)*cos(q4)*cos(q1 + pi/4)), 
                       0], 
                      [L1*cos(q1 + pi/4) + L2*cos(q2)*cos(q1 + pi/4) - L3*sin(q2)*cos(q3)*cos(q1 + pi/4) - L3*sin(q3)*sin(q1 + pi/4) - L4*sin(q2)*sin(q4)*cos(q3)*cos(q1 + pi/4) - L4*sin(q3)*sin(q4)*sin(q1 + pi/4) + L4*cos(q2)*cos(q4)*cos(q1 + pi/4) + L5*sin(q2)*sin(q3)*sin(q5)*cos(q1 + pi/4) - L5*sin(q2)*cos(q3)*cos(q4)*cos(q5)*cos(q1 + pi/4) - L5*sin(q3)*sin(q1 + pi/4)*cos(q4)*cos(q5) - L5*sin(q4)*cos(q2)*cos(q5)*cos(q1 + pi/4) - L5*sin(q5)*sin(q1 + pi/4)*cos(q3) + L6*sin(q2)*sin(q3)*sin(q5)*sin(q6)*cos(q1 + pi/4) - L6*sin(q2)*sin(q4)*cos(q3)*cos(q6)*cos(q1 + pi/4) - L6*sin(q2)*sin(q6)*cos(q3)*cos(q4)*cos(q5)*cos(q1 + pi/4) - L6*sin(q3)*sin(q4)*sin(q1 + pi/4)*cos(q6) - L6*sin(q3)*sin(q6)*sin(q1 + pi/4)*cos(q4)*cos(q5) - L6*sin(q4)*sin(q6)*cos(q2)*cos(q5)*cos(q1 + pi/4) - L6*sin(q5)*sin(q6)*sin(q1 + pi/4)*cos(q3) + L6*cos(q2)*cos(q4)*cos(q6)*cos(q1 + pi/4), 
                       (-L2*sin(q2) - L3*cos(q2)*cos(q3) - L4*sin(q2)*cos(q4) - L4*sin(q4)*cos(q2)*cos(q3) + L5*sin(q2)*sin(q4)*cos(q5) + L5*sin(q3)*sin(q5)*cos(q2) - L5*cos(q2)*cos(q3)*cos(q4)*cos(q5) + L6*sin(q2)*sin(q4)*sin(q6)*cos(q5) - L6*sin(q2)*cos(q4)*cos(q6) + L6*sin(q3)*sin(q5)*sin(q6)*cos(q2) - L6*sin(q4)*cos(q2)*cos(q3)*cos(q6) - L6*sin(q6)*cos(q2)*cos(q3)*cos(q4)*cos(q5))*sin(q1 + pi/4), 
                       L3*sin(q2)*sin(q3)*sin(q1 + pi/4) + L3*cos(q3)*cos(q1 + pi/4) + L4*sin(q2)*sin(q3)*sin(q4)*sin(q1 + pi/4) + L4*sin(q4)*cos(q3)*cos(q1 + pi/4) + L5*sin(q2)*sin(q3)*sin(q1 + pi/4)*cos(q4)*cos(q5) + L5*sin(q2)*sin(q5)*sin(q1 + pi/4)*cos(q3) - L5*sin(q3)*sin(q5)*cos(q1 + pi/4) + L5*cos(q3)*cos(q4)*cos(q5)*cos(q1 + pi/4) + L6*sin(q2)*sin(q3)*sin(q4)*sin(q1 + pi/4)*cos(q6) + L6*sin(q2)*sin(q3)*sin(q6)*sin(q1 + pi/4)*cos(q4)*cos(q5) + L6*sin(q2)*sin(q5)*sin(q6)*sin(q1 + pi/4)*cos(q3) - L6*sin(q3)*sin(q5)*sin(q6)*cos(q1 + pi/4) + L6*sin(q4)*cos(q3)*cos(q6)*cos(q1 + pi/4) + L6*sin(q6)*cos(q3)*cos(q4)*cos(q5)*cos(q1 + pi/4), 
                       -L4*sin(q2)*sin(q1 + pi/4)*cos(q3)*cos(q4) + L4*sin(q3)*cos(q4)*cos(q1 + pi/4) - L4*sin(q4)*sin(q1 + pi/4)*cos(q2) + L5*sin(q2)*sin(q4)*sin(q1 + pi/4)*cos(q3)*cos(q5) - L5*sin(q3)*sin(q4)*cos(q5)*cos(q1 + pi/4) - L5*sin(q1 + pi/4)*cos(q2)*cos(q4)*cos(q5) + L6*sin(q2)*sin(q4)*sin(q6)*sin(q1 + pi/4)*cos(q3)*cos(q5) - L6*sin(q2)*sin(q1 + pi/4)*cos(q3)*cos(q4)*cos(q6) - L6*sin(q3)*sin(q4)*sin(q6)*cos(q5)*cos(q1 + pi/4) + L6*sin(q3)*cos(q4)*cos(q6)*cos(q1 + pi/4) - L6*sin(q4)*sin(q1 + pi/4)*cos(q2)*cos(q6) - L6*sin(q6)*sin(q1 + pi/4)*cos(q2)*cos(q4)*cos(q5), 
                       L5*sin(q2)*sin(q3)*sin(q1 + pi/4)*cos(q5) + L5*sin(q2)*sin(q5)*sin(q1 + pi/4)*cos(q3)*cos(q4) - L5*sin(q3)*sin(q5)*cos(q4)*cos(q1 + pi/4) + L5*sin(q4)*sin(q5)*sin(q1 + pi/4)*cos(q2) + L5*cos(q3)*cos(q5)*cos(q1 + pi/4) + L6*sin(q2)*sin(q3)*sin(q6)*sin(q1 + pi/4)*cos(q5) + L6*sin(q2)*sin(q5)*sin(q6)*sin(q1 + pi/4)*cos(q3)*cos(q4) - L6*sin(q3)*sin(q5)*sin(q6)*cos(q4)*cos(q1 + pi/4) + L6*sin(q4)*sin(q5)*sin(q6)*sin(q1 + pi/4)*cos(q2) + L6*sin(q6)*cos(q3)*cos(q5)*cos(q1 + pi/4), 
                       L6*(sin(q2)*sin(q3)*sin(q5)*sin(q1 + pi/4)*cos(q6) + sin(q2)*sin(q4)*sin(q6)*sin(q1 + pi/4)*cos(q3) - sin(q2)*sin(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5)*cos(q6) - sin(q3)*sin(q4)*sin(q6)*cos(q1 + pi/4) + sin(q3)*cos(q4)*cos(q5)*cos(q6)*cos(q1 + pi/4) - sin(q4)*sin(q1 + pi/4)*cos(q2)*cos(q5)*cos(q6) + sin(q5)*cos(q3)*cos(q6)*cos(q1 + pi/4) - sin(q6)*sin(q1 + pi/4)*cos(q2)*cos(q4)), 
                       0], 
                      [0, 
                       -L2*cos(q2) + L3*sin(q2)*cos(q3) + L4*(sin(q2)*sin(q4)*cos(q3) - cos(q2)*cos(q4)) + L5*((sin(q2)*cos(q3)*cos(q4) + sin(q4)*cos(q2))*cos(q5) - sin(q2)*sin(q3)*sin(q5)) + L6*(((sin(q2)*cos(q3)*cos(q4) + sin(q4)*cos(q2))*cos(q5) - sin(q2)*sin(q3)*sin(q5))*sin(q6) + (sin(q2)*sin(q4)*cos(q3) - cos(q2)*cos(q4))*cos(q6)), 
                       (L3*sin(q3) + L4*sin(q3)*sin(q4) + L5*(sin(q3)*cos(q4)*cos(q5) + sin(q5)*cos(q3)) + L6*((sin(q3)*cos(q4)*cos(q5) + sin(q5)*cos(q3))*sin(q6) + sin(q3)*sin(q4)*cos(q6)))*cos(q2), 
                       L4*(sin(q2)*sin(q4) - cos(q2)*cos(q3)*cos(q4)) + L5*(sin(q2)*cos(q4) + sin(q4)*cos(q2)*cos(q3))*cos(q5) + L6*((sin(q2)*sin(q4) - cos(q2)*cos(q3)*cos(q4))*cos(q6) + (sin(q2)*cos(q4) + sin(q4)*cos(q2)*cos(q3))*sin(q6)*cos(q5)), 
                       -(L5 + L6*sin(q6))*((sin(q2)*sin(q4) - cos(q2)*cos(q3)*cos(q4))*sin(q5) - sin(q3)*cos(q2)*cos(q5)), 
                       L6*(((sin(q2)*sin(q4) - cos(q2)*cos(q3)*cos(q4))*cos(q5) + sin(q3)*sin(q5)*cos(q2))*cos(q6) + (sin(q2)*cos(q4) + sin(q4)*cos(q2)*cos(q3))*sin(q6)), 
                       0]])     
        return z
    
    def jacobi_all(self, q):
        """ヤコビ行列を全部計算"""
        jWo = self.jacobi_Wo(q)
        jBL = self.jacobi_BL(q)
        j0 = self.jacobi_0(q)
        j1 = self.jacobi_1(q)
        j2 = self.jacobi_2(q)
        j3 = self.jacobi_3(q)
        j4 = self.jacobi_4(q)
        j5 = self.jacobi_5(q)
        j6 = self.jacobi_6(q)
        j7 = self.jacobi_7(q)
        jGL = self.jacobi_GL(q)
        return [jWo, jBL, j0, j1, j2, j3, j4, j5, j6, j7, jGL]
    
    
    def djacobi_Wo(self, q, dq):
        L, h, H, L0, L1, L2, L3, L4, L5, L6 = self.L, self.h, self.H, self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6
        q1, q2, q3, q4, q5, q6, q7 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0], q[6, 0]
        dq1, dq2, dq3, dq4, dq5, dq6, dq7 = dq[0, 0], dq[1, 0], dq[2, 0], dq[3, 0], dq[4, 0], dq[5, 0], dq[6, 0]
        z = np.array([[0, 
                       0, 
                       0, 
                       0, 
                       0, 
                       0, 
                       0], 
                      [0, 
                       0, 
                       0, 
                       0, 
                       0, 
                       0, 
                       0], 
                      [0, 
                       0, 
                       0, 
                       0, 
                       0, 
                       0, 
                       0]])
        return z
    
    def djacobi_BL(self, q, dq):
        L, h, H, L0, L1, L2, L3, L4, L5, L6 = self.L, self.h, self.H, self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6
        q1, q2, q3, q4, q5, q6, q7 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0], q[6, 0]
        dq1, dq2, dq3, dq4, dq5, dq6, dq7 = dq[0, 0], dq[1, 0], dq[2, 0], dq[3, 0], dq[4, 0], dq[5, 0], dq[6, 0]
        z = np.array([[0, 
                       0, 
                       0, 
                       0, 
                       0, 
                       0, 
                       0], 
                      [0, 
                       0, 
                       0, 
                       0, 
                       0, 
                       0, 
                       0], 
                      [0, 
                       0, 
                       0, 
                       0, 
                       0, 
                       0, 
                       0]])
        return z
    
    def djacobi_0(self, q, dq):
        L, h, H, L0, L1, L2, L3, L4, L5, L6 = self.L, self.h, self.H, self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6
        q1, q2, q3, q4, q5, q6, q7 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0], q[6, 0]
        dq1, dq2, dq3, dq4, dq5, dq6, dq7 = dq[0, 0], dq[1, 0], dq[2, 0], dq[3, 0], dq[4, 0], dq[5, 0], dq[6, 0]
        z = np.array([[0, 0, 0, 0, 0, 0, 0], 
                      [0, 0, 0, 0, 0, 0, 0], 
                      [0, 0, 0, 0, 0, 0, 0]])
        return z
    
    def djacobi_1(self, q, dq):
        L, h, H, L0, L1, L2, L3, L4, L5, L6 = self.L, self.h, self.H, self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6
        q1, q2, q3, q4, q5, q6, q7 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0], q[6, 0]
        dq1, dq2, dq3, dq4, dq5, dq6, dq7 = dq[0, 0], dq[1, 0], dq[2, 0], dq[3, 0], dq[4, 0], dq[5, 0], dq[6, 0]
        z = np.array([[0, 0, 0, 0, 0, 0, 0], 
                      [0, 0, 0, 0, 0, 0, 0], 
                      [0, 0, 0, 0, 0, 0, 0]])
        return z
    
    def djacobi_2(self, q, dq):
        L, h, H, L0, L1, L2, L3, L4, L5, L6 = self.L, self.h, self.H, self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6
        q1, q2, q3, q4, q5, q6, q7 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0], q[6, 0]
        dq1, dq2, dq3, dq4, dq5, dq6, dq7 = dq[0, 0], dq[1, 0], dq[2, 0], dq[3, 0], dq[4, 0], dq[5, 0], dq[6, 0]
        z = np.array([[-L1*dq1*cos(q1 + pi/4), 
                       0, 
                       0, 
                       0, 
                       0, 
                       0, 
                       0], 
                      [-L1*dq1*sin(q1 + pi/4), 
                       0, 
                       0, 
                       0, 
                       0, 
                       0, 
                       0], 
                      [0, 
                       0, 
                       0, 
                       0, 
                       0, 
                       0, 
                       0]])
        return z
    
    def djacobi_3(self, q, dq):
        L, h, H, L0, L1, L2, L3, L4, L5, L6 = self.L, self.h, self.H, self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6
        q1, q2, q3, q4, q5, q6, q7 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0], q[6, 0]
        dq1, dq2, dq3, dq4, dq5, dq6, dq7 = dq[0, 0], dq[1, 0], dq[2, 0], dq[3, 0], dq[4, 0], dq[5, 0], dq[6, 0]
        z = np.array([[L2*dq2*sin(q1 + pi/4)*sin(q2) + dq1*(-L1 - L2*cos(q2))*cos(q1 + pi/4), 
                       L2*dq1*sin(q1 + pi/4)*sin(q2) - L2*dq2*cos(q1 + pi/4)*cos(q2), 
                       0, 
                       0, 
                       0, 
                       0, 
                       0], 
                      [-L2*dq2*sin(q2)*cos(q1 + pi/4) - dq1*(L1 + L2*cos(q2))*sin(q1 + pi/4), 
                       -L2*dq1*sin(q2)*cos(q1 + pi/4) - L2*dq2*sin(q1 + pi/4)*cos(q2), 
                       0, 
                       0, 
                       0, 
                       0, 
                       0], 
                      [0, 
                       L2*dq2*sin(q2), 
                       0, 
                       0, 
                       0, 
                       0, 
                       0]])
        return z
    
    def djacobi_4(self, q, dq):
        L, h, H, L0, L1, L2, L3, L4, L5, L6 = self.L, self.h, self.H, self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6
        q1, q2, q3, q4, q5, q6, q7 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0], q[6, 0]
        dq1, dq2, dq3, dq4, dq5, dq6, dq7 = dq[0, 0], dq[1, 0], dq[2, 0], dq[3, 0], dq[4, 0], dq[5, 0], dq[6, 0]
        z = np.array([[-L1*dq1*cos(q1 + pi/4) - L2*dq1*cos(q1 + pi/4)*cos(q2) + L2*dq2*sin(q1 + pi/4)*sin(q2) + L3*dq1*sin(q1 + pi/4)*sin(q3) + L3*dq1*sin(q2)*cos(q1 + pi/4)*cos(q3) + L3*dq2*sin(q1 + pi/4)*cos(q2)*cos(q3) - L3*dq3*sin(q1 + pi/4)*sin(q2)*sin(q3) - L3*dq3*cos(q1 + pi/4)*cos(q3), 
                       -dq1*(-L2*sin(q2) - L3*cos(q2)*cos(q3))*sin(q1 + pi/4) + (-L2*dq2*cos(q2) + L3*dq2*sin(q2)*cos(q3) + L3*dq3*sin(q3)*cos(q2))*cos(q1 + pi/4), 
                       L3*(-dq1*sin(q1 + pi/4)*sin(q2)*sin(q3) - dq1*cos(q1 + pi/4)*cos(q3) + dq2*sin(q3)*cos(q1 + pi/4)*cos(q2) + dq3*sin(q1 + pi/4)*sin(q3) + dq3*sin(q2)*cos(q1 + pi/4)*cos(q3)), 
                       0, 
                       0, 
                       0, 
                       0], 
                      [-L1*dq1*sin(q1 + pi/4) - L2*dq1*sin(q1 + pi/4)*cos(q2) - L2*dq2*sin(q2)*cos(q1 + pi/4) + L3*dq1*sin(q1 + pi/4)*sin(q2)*cos(q3) - L3*dq1*sin(q3)*cos(q1 + pi/4) - L3*dq2*cos(q1 + pi/4)*cos(q2)*cos(q3) - L3*dq3*sin(q1 + pi/4)*cos(q3) + L3*dq3*sin(q2)*sin(q3)*cos(q1 + pi/4), 
                       dq1*(-L2*sin(q2) - L3*cos(q2)*cos(q3))*cos(q1 + pi/4) + (-L2*dq2*cos(q2) + L3*dq2*sin(q2)*cos(q3) + L3*dq3*sin(q3)*cos(q2))*sin(q1 + pi/4), 
                       L3*(-dq1*sin(q1 + pi/4)*cos(q3) + dq1*sin(q2)*sin(q3)*cos(q1 + pi/4) + dq2*sin(q1 + pi/4)*sin(q3)*cos(q2) + dq3*sin(q1 + pi/4)*sin(q2)*cos(q3) - dq3*sin(q3)*cos(q1 + pi/4)), 
                       0, 
                       0, 
                       0, 
                       0], 
                      [0, 
                       L2*dq2*sin(q2) + L3*dq2*cos(q2)*cos(q3) - L3*dq3*sin(q2)*sin(q3), 
                       -L3*dq2*sin(q2)*sin(q3) + L3*dq3*cos(q2)*cos(q3), 
                       0, 
                       0, 
                       0, 
                       0]])
        return z
    
    def djacobi_5(self, q, dq):
        L, h, H, L0, L1, L2, L3, L4, L5, L6 = self.L, self.h, self.H, self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6
        q1, q2, q3, q4, q5, q6, q7 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0], q[6, 0]
        dq1, dq2, dq3, dq4, dq5, dq6, dq7 = dq[0, 0], dq[1, 0], dq[2, 0], dq[3, 0], dq[4, 0], dq[5, 0], dq[6, 0]
        z = np.array([[-L1*dq1*cos(q1 + pi/4) - L2*dq1*cos(q1 + pi/4)*cos(q2) + L2*dq2*sin(q1 + pi/4)*sin(q2) + L3*dq1*sin(q1 + pi/4)*sin(q3) + L3*dq1*sin(q2)*cos(q1 + pi/4)*cos(q3) + L3*dq2*sin(q1 + pi/4)*cos(q2)*cos(q3) - L3*dq3*sin(q1 + pi/4)*sin(q2)*sin(q3) - L3*dq3*cos(q1 + pi/4)*cos(q3) + L4*dq1*sin(q1 + pi/4)*sin(q3)*sin(q4) + L4*dq1*sin(q2)*sin(q4)*cos(q1 + pi/4)*cos(q3) - L4*dq1*cos(q1 + pi/4)*cos(q2)*cos(q4) + L4*dq2*sin(q1 + pi/4)*sin(q2)*cos(q4) + L4*dq2*sin(q1 + pi/4)*sin(q4)*cos(q2)*cos(q3) - L4*dq3*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q4) - L4*dq3*sin(q4)*cos(q1 + pi/4)*cos(q3) + L4*dq4*sin(q1 + pi/4)*sin(q2)*cos(q3)*cos(q4) + L4*dq4*sin(q1 + pi/4)*sin(q4)*cos(q2) - L4*dq4*sin(q3)*cos(q1 + pi/4)*cos(q4), 
                       -dq1*(-L2*sin(q2) - L3*cos(q2)*cos(q3) - L4*sin(q2)*cos(q4) - L4*sin(q4)*cos(q2)*cos(q3))*sin(q1 + pi/4) + (-L2*dq2*cos(q2) + L3*dq2*sin(q2)*cos(q3) + L3*dq3*sin(q3)*cos(q2) + L4*dq2*sin(q2)*sin(q4)*cos(q3) - L4*dq2*cos(q2)*cos(q4) + L4*dq3*sin(q3)*sin(q4)*cos(q2) + L4*dq4*sin(q2)*sin(q4) - L4*dq4*cos(q2)*cos(q3)*cos(q4))*cos(q1 + pi/4), 
                       -L3*dq1*sin(q1 + pi/4)*sin(q2)*sin(q3) - L3*dq1*cos(q1 + pi/4)*cos(q3) + L3*dq2*sin(q3)*cos(q1 + pi/4)*cos(q2) + L3*dq3*sin(q1 + pi/4)*sin(q3) + L3*dq3*sin(q2)*cos(q1 + pi/4)*cos(q3) - L4*dq1*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q4) - L4*dq1*sin(q4)*cos(q1 + pi/4)*cos(q3) + L4*dq2*sin(q3)*sin(q4)*cos(q1 + pi/4)*cos(q2) + L4*dq3*sin(q1 + pi/4)*sin(q3)*sin(q4) + L4*dq3*sin(q2)*sin(q4)*cos(q1 + pi/4)*cos(q3) - L4*dq4*sin(q1 + pi/4)*cos(q3)*cos(q4) + L4*dq4*sin(q2)*sin(q3)*cos(q1 + pi/4)*cos(q4), 
                       -L4*(-dq1*sin(q1 + pi/4)*sin(q2)*cos(q3)*cos(q4) - dq1*sin(q1 + pi/4)*sin(q4)*cos(q2) + dq1*sin(q3)*cos(q1 + pi/4)*cos(q4) - dq2*sin(q2)*sin(q4)*cos(q1 + pi/4) + dq2*cos(q1 + pi/4)*cos(q2)*cos(q3)*cos(q4) + dq3*sin(q1 + pi/4)*cos(q3)*cos(q4) - dq3*sin(q2)*sin(q3)*cos(q1 + pi/4)*cos(q4) - dq4*sin(q1 + pi/4)*sin(q3)*sin(q4) - dq4*sin(q2)*sin(q4)*cos(q1 + pi/4)*cos(q3) + dq4*cos(q1 + pi/4)*cos(q2)*cos(q4)), 
                       0, 
                       0, 
                       0], 
                      [-L1*dq1*sin(q1 + pi/4) - L2*dq1*sin(q1 + pi/4)*cos(q2) - L2*dq2*sin(q2)*cos(q1 + pi/4) + L3*dq1*sin(q1 + pi/4)*sin(q2)*cos(q3) - L3*dq1*sin(q3)*cos(q1 + pi/4) - L3*dq2*cos(q1 + pi/4)*cos(q2)*cos(q3) - L3*dq3*sin(q1 + pi/4)*cos(q3) + L3*dq3*sin(q2)*sin(q3)*cos(q1 + pi/4) + L4*dq1*sin(q1 + pi/4)*sin(q2)*sin(q4)*cos(q3) - L4*dq1*sin(q1 + pi/4)*cos(q2)*cos(q4) - L4*dq1*sin(q3)*sin(q4)*cos(q1 + pi/4) - L4*dq2*sin(q2)*cos(q1 + pi/4)*cos(q4) - L4*dq2*sin(q4)*cos(q1 + pi/4)*cos(q2)*cos(q3) - L4*dq3*sin(q1 + pi/4)*sin(q4)*cos(q3) + L4*dq3*sin(q2)*sin(q3)*sin(q4)*cos(q1 + pi/4) - L4*dq4*sin(q1 + pi/4)*sin(q3)*cos(q4) - L4*dq4*sin(q2)*cos(q1 + pi/4)*cos(q3)*cos(q4) - L4*dq4*sin(q4)*cos(q1 + pi/4)*cos(q2), 
                       dq1*(-L2*sin(q2) - L3*cos(q2)*cos(q3) - L4*sin(q2)*cos(q4) - L4*sin(q4)*cos(q2)*cos(q3))*cos(q1 + pi/4) + (-L2*dq2*cos(q2) + L3*dq2*sin(q2)*cos(q3) + L3*dq3*sin(q3)*cos(q2) + L4*dq2*sin(q2)*sin(q4)*cos(q3) - L4*dq2*cos(q2)*cos(q4) + L4*dq3*sin(q3)*sin(q4)*cos(q2) + L4*dq4*sin(q2)*sin(q4) - L4*dq4*cos(q2)*cos(q3)*cos(q4))*sin(q1 + pi/4), 
                       -L3*dq1*sin(q1 + pi/4)*cos(q3) + L3*dq1*sin(q2)*sin(q3)*cos(q1 + pi/4) + L3*dq2*sin(q1 + pi/4)*sin(q3)*cos(q2) + L3*dq3*sin(q1 + pi/4)*sin(q2)*cos(q3) - L3*dq3*sin(q3)*cos(q1 + pi/4) - L4*dq1*sin(q1 + pi/4)*sin(q4)*cos(q3) + L4*dq1*sin(q2)*sin(q3)*sin(q4)*cos(q1 + pi/4) + L4*dq2*sin(q1 + pi/4)*sin(q3)*sin(q4)*cos(q2) + L4*dq3*sin(q1 + pi/4)*sin(q2)*sin(q4)*cos(q3) - L4*dq3*sin(q3)*sin(q4)*cos(q1 + pi/4) + L4*dq4*sin(q1 + pi/4)*sin(q2)*sin(q3)*cos(q4) + L4*dq4*cos(q1 + pi/4)*cos(q3)*cos(q4), 
                       L4*(-dq1*sin(q1 + pi/4)*sin(q3)*cos(q4) - dq1*sin(q2)*cos(q1 + pi/4)*cos(q3)*cos(q4) - dq1*sin(q4)*cos(q1 + pi/4)*cos(q2) + dq2*sin(q1 + pi/4)*sin(q2)*sin(q4) - dq2*sin(q1 + pi/4)*cos(q2)*cos(q3)*cos(q4) + dq3*sin(q1 + pi/4)*sin(q2)*sin(q3)*cos(q4) + dq3*cos(q1 + pi/4)*cos(q3)*cos(q4) + dq4*sin(q1 + pi/4)*sin(q2)*sin(q4)*cos(q3) - dq4*sin(q1 + pi/4)*cos(q2)*cos(q4) - dq4*sin(q3)*sin(q4)*cos(q1 + pi/4)), 
                       0, 
                       0, 
                       0], 
                      [0, 
                       L2*dq2*sin(q2) + L3*dq2*cos(q2)*cos(q3) - L3*dq3*sin(q2)*sin(q3) + L4*(dq2*sin(q2)*cos(q4) + dq2*sin(q4)*cos(q2)*cos(q3) - dq3*sin(q2)*sin(q3)*sin(q4) + dq4*sin(q2)*cos(q3)*cos(q4) + dq4*sin(q4)*cos(q2)), 
                       L4*dq4*sin(q3)*cos(q2)*cos(q4) - dq2*(L3 + L4*sin(q4))*sin(q2)*sin(q3) + dq3*(L3 + L4*sin(q4))*cos(q2)*cos(q3), 
                       L4*(dq2*sin(q2)*cos(q3)*cos(q4) + dq2*sin(q4)*cos(q2) + dq3*sin(q3)*cos(q2)*cos(q4) + dq4*sin(q2)*cos(q4) + dq4*sin(q4)*cos(q2)*cos(q3)), 
                       0, 
                       0, 
                       0]])
        return z
    
    def djacobi_6(self, q, dq):
        L, h, H, L0, L1, L2, L3, L4, L5, L6 = self.L, self.h, self.H, self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6
        q1, q2, q3, q4, q5, q6, q7 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0], q[6, 0]
        dq1, dq2, dq3, dq4, dq5, dq6, dq7 = dq[0, 0], dq[1, 0], dq[2, 0], dq[3, 0], dq[4, 0], dq[5, 0], dq[6, 0]
        z = np.array([[-L1*dq1*cos(q1 + pi/4) - L2*dq1*cos(q1 + pi/4)*cos(q2) + L2*dq2*sin(q1 + pi/4)*sin(q2) + L3*dq1*sin(q1 + pi/4)*sin(q3) + L3*dq1*sin(q2)*cos(q1 + pi/4)*cos(q3) + L3*dq2*sin(q1 + pi/4)*cos(q2)*cos(q3) - L3*dq3*sin(q1 + pi/4)*sin(q2)*sin(q3) - L3*dq3*cos(q1 + pi/4)*cos(q3) + L4*dq1*sin(q1 + pi/4)*sin(q3)*sin(q4) + L4*dq1*sin(q2)*sin(q4)*cos(q1 + pi/4)*cos(q3) - L4*dq1*cos(q1 + pi/4)*cos(q2)*cos(q4) + L4*dq2*sin(q1 + pi/4)*sin(q2)*cos(q4) + L4*dq2*sin(q1 + pi/4)*sin(q4)*cos(q2)*cos(q3) - L4*dq3*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q4) - L4*dq3*sin(q4)*cos(q1 + pi/4)*cos(q3) + L4*dq4*sin(q1 + pi/4)*sin(q2)*cos(q3)*cos(q4) + L4*dq4*sin(q1 + pi/4)*sin(q4)*cos(q2) - L4*dq4*sin(q3)*cos(q1 + pi/4)*cos(q4) + L5*dq1*sin(q1 + pi/4)*sin(q3)*cos(q4)*cos(q5) + L5*dq1*sin(q1 + pi/4)*sin(q5)*cos(q3) - L5*dq1*sin(q2)*sin(q3)*sin(q5)*cos(q1 + pi/4) + L5*dq1*sin(q2)*cos(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) + L5*dq1*sin(q4)*cos(q1 + pi/4)*cos(q2)*cos(q5) - L5*dq2*sin(q1 + pi/4)*sin(q2)*sin(q4)*cos(q5) - L5*dq2*sin(q1 + pi/4)*sin(q3)*sin(q5)*cos(q2) + L5*dq2*sin(q1 + pi/4)*cos(q2)*cos(q3)*cos(q4)*cos(q5) - L5*dq3*sin(q1 + pi/4)*sin(q2)*sin(q3)*cos(q4)*cos(q5) - L5*dq3*sin(q1 + pi/4)*sin(q2)*sin(q5)*cos(q3) + L5*dq3*sin(q3)*sin(q5)*cos(q1 + pi/4) - L5*dq3*cos(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) - L5*dq4*sin(q1 + pi/4)*sin(q2)*sin(q4)*cos(q3)*cos(q5) + L5*dq4*sin(q1 + pi/4)*cos(q2)*cos(q4)*cos(q5) + L5*dq4*sin(q3)*sin(q4)*cos(q1 + pi/4)*cos(q5) - L5*dq5*sin(q1 + pi/4)*sin(q2)*sin(q3)*cos(q5) - L5*dq5*sin(q1 + pi/4)*sin(q2)*sin(q5)*cos(q3)*cos(q4) - L5*dq5*sin(q1 + pi/4)*sin(q4)*sin(q5)*cos(q2) + L5*dq5*sin(q3)*sin(q5)*cos(q1 + pi/4)*cos(q4) - L5*dq5*cos(q1 + pi/4)*cos(q3)*cos(q5), 
                       -dq1*(-L2*sin(q2) - L3*cos(q2)*cos(q3) - L4*sin(q2)*cos(q4) - L4*sin(q4)*cos(q2)*cos(q3) + L5*sin(q2)*sin(q4)*cos(q5) + L5*sin(q3)*sin(q5)*cos(q2) - L5*cos(q2)*cos(q3)*cos(q4)*cos(q5))*sin(q1 + pi/4) + (-L2*dq2*cos(q2) + L3*dq2*sin(q2)*cos(q3) + L3*dq3*sin(q3)*cos(q2) + L4*dq2*sin(q2)*sin(q4)*cos(q3) - L4*dq2*cos(q2)*cos(q4) + L4*dq3*sin(q3)*sin(q4)*cos(q2) + L4*dq4*sin(q2)*sin(q4) - L4*dq4*cos(q2)*cos(q3)*cos(q4) - L5*dq2*sin(q2)*sin(q3)*sin(q5) + L5*dq2*sin(q2)*cos(q3)*cos(q4)*cos(q5) + L5*dq2*sin(q4)*cos(q2)*cos(q5) + L5*dq3*sin(q3)*cos(q2)*cos(q4)*cos(q5) + L5*dq3*sin(q5)*cos(q2)*cos(q3) + L5*dq4*sin(q2)*cos(q4)*cos(q5) + L5*dq4*sin(q4)*cos(q2)*cos(q3)*cos(q5) - L5*dq5*sin(q2)*sin(q4)*sin(q5) + L5*dq5*sin(q3)*cos(q2)*cos(q5) + L5*dq5*sin(q5)*cos(q2)*cos(q3)*cos(q4))*cos(q1 + pi/4), 
                       -L3*dq1*sin(q1 + pi/4)*sin(q2)*sin(q3) - L3*dq1*cos(q1 + pi/4)*cos(q3) + L3*dq2*sin(q3)*cos(q1 + pi/4)*cos(q2) + L3*dq3*sin(q1 + pi/4)*sin(q3) + L3*dq3*sin(q2)*cos(q1 + pi/4)*cos(q3) - L4*dq1*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q4) - L4*dq1*sin(q4)*cos(q1 + pi/4)*cos(q3) + L4*dq2*sin(q3)*sin(q4)*cos(q1 + pi/4)*cos(q2) + L4*dq3*sin(q1 + pi/4)*sin(q3)*sin(q4) + L4*dq3*sin(q2)*sin(q4)*cos(q1 + pi/4)*cos(q3) - L4*dq4*sin(q1 + pi/4)*cos(q3)*cos(q4) + L4*dq4*sin(q2)*sin(q3)*cos(q1 + pi/4)*cos(q4) - L5*dq1*sin(q1 + pi/4)*sin(q2)*sin(q3)*cos(q4)*cos(q5) - L5*dq1*sin(q1 + pi/4)*sin(q2)*sin(q5)*cos(q3) + L5*dq1*sin(q3)*sin(q5)*cos(q1 + pi/4) - L5*dq1*cos(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) + L5*dq2*sin(q3)*cos(q1 + pi/4)*cos(q2)*cos(q4)*cos(q5) + L5*dq2*sin(q5)*cos(q1 + pi/4)*cos(q2)*cos(q3) + L5*dq3*sin(q1 + pi/4)*sin(q3)*cos(q4)*cos(q5) + L5*dq3*sin(q1 + pi/4)*sin(q5)*cos(q3) - L5*dq3*sin(q2)*sin(q3)*sin(q5)*cos(q1 + pi/4) + L5*dq3*sin(q2)*cos(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) + L5*dq4*sin(q1 + pi/4)*sin(q4)*cos(q3)*cos(q5) - L5*dq4*sin(q2)*sin(q3)*sin(q4)*cos(q1 + pi/4)*cos(q5) + L5*dq5*sin(q1 + pi/4)*sin(q3)*cos(q5) + L5*dq5*sin(q1 + pi/4)*sin(q5)*cos(q3)*cos(q4) - L5*dq5*sin(q2)*sin(q3)*sin(q5)*cos(q1 + pi/4)*cos(q4) + L5*dq5*sin(q2)*cos(q1 + pi/4)*cos(q3)*cos(q5), 
                       L4*dq1*sin(q1 + pi/4)*sin(q2)*cos(q3)*cos(q4) + L4*dq1*sin(q1 + pi/4)*sin(q4)*cos(q2) - L4*dq1*sin(q3)*cos(q1 + pi/4)*cos(q4) + L4*dq2*sin(q2)*sin(q4)*cos(q1 + pi/4) - L4*dq2*cos(q1 + pi/4)*cos(q2)*cos(q3)*cos(q4) - L4*dq3*sin(q1 + pi/4)*cos(q3)*cos(q4) + L4*dq3*sin(q2)*sin(q3)*cos(q1 + pi/4)*cos(q4) + L4*dq4*sin(q1 + pi/4)*sin(q3)*sin(q4) + L4*dq4*sin(q2)*sin(q4)*cos(q1 + pi/4)*cos(q3) - L4*dq4*cos(q1 + pi/4)*cos(q2)*cos(q4) - L5*dq1*sin(q1 + pi/4)*sin(q2)*sin(q4)*cos(q3)*cos(q5) + L5*dq1*sin(q1 + pi/4)*cos(q2)*cos(q4)*cos(q5) + L5*dq1*sin(q3)*sin(q4)*cos(q1 + pi/4)*cos(q5) + L5*dq2*sin(q2)*cos(q1 + pi/4)*cos(q4)*cos(q5) + L5*dq2*sin(q4)*cos(q1 + pi/4)*cos(q2)*cos(q3)*cos(q5) + L5*dq3*sin(q1 + pi/4)*sin(q4)*cos(q3)*cos(q5) - L5*dq3*sin(q2)*sin(q3)*sin(q4)*cos(q1 + pi/4)*cos(q5) + L5*dq4*sin(q1 + pi/4)*sin(q3)*cos(q4)*cos(q5) + L5*dq4*sin(q2)*cos(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) + L5*dq4*sin(q4)*cos(q1 + pi/4)*cos(q2)*cos(q5) - L5*dq5*sin(q1 + pi/4)*sin(q3)*sin(q4)*sin(q5) - L5*dq5*sin(q2)*sin(q4)*sin(q5)*cos(q1 + pi/4)*cos(q3) + L5*dq5*sin(q5)*cos(q1 + pi/4)*cos(q2)*cos(q4), 
                       L5*(-dq1*sin(q1 + pi/4)*sin(q2)*sin(q3)*cos(q5) - dq1*sin(q1 + pi/4)*sin(q2)*sin(q5)*cos(q3)*cos(q4) - dq1*sin(q1 + pi/4)*sin(q4)*sin(q5)*cos(q2) + dq1*sin(q3)*sin(q5)*cos(q1 + pi/4)*cos(q4) - dq1*cos(q1 + pi/4)*cos(q3)*cos(q5) - dq2*sin(q2)*sin(q4)*sin(q5)*cos(q1 + pi/4) + dq2*sin(q3)*cos(q1 + pi/4)*cos(q2)*cos(q5) + dq2*sin(q5)*cos(q1 + pi/4)*cos(q2)*cos(q3)*cos(q4) + dq3*sin(q1 + pi/4)*sin(q3)*cos(q5) + dq3*sin(q1 + pi/4)*sin(q5)*cos(q3)*cos(q4) - dq3*sin(q2)*sin(q3)*sin(q5)*cos(q1 + pi/4)*cos(q4) + dq3*sin(q2)*cos(q1 + pi/4)*cos(q3)*cos(q5) - dq4*sin(q1 + pi/4)*sin(q3)*sin(q4)*sin(q5) - dq4*sin(q2)*sin(q4)*sin(q5)*cos(q1 + pi/4)*cos(q3) + dq4*sin(q5)*cos(q1 + pi/4)*cos(q2)*cos(q4) + dq5*sin(q1 + pi/4)*sin(q3)*cos(q4)*cos(q5) + dq5*sin(q1 + pi/4)*sin(q5)*cos(q3) - dq5*sin(q2)*sin(q3)*sin(q5)*cos(q1 + pi/4) + dq5*sin(q2)*cos(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) + dq5*sin(q4)*cos(q1 + pi/4)*cos(q2)*cos(q5)), 
                       0, 
                       0], 
                      [-L1*dq1*sin(q1 + pi/4) - L2*dq1*sin(q1 + pi/4)*cos(q2) - L2*dq2*sin(q2)*cos(q1 + pi/4) + L3*dq1*sin(q1 + pi/4)*sin(q2)*cos(q3) - L3*dq1*sin(q3)*cos(q1 + pi/4) - L3*dq2*cos(q1 + pi/4)*cos(q2)*cos(q3) - L3*dq3*sin(q1 + pi/4)*cos(q3) + L3*dq3*sin(q2)*sin(q3)*cos(q1 + pi/4) + L4*dq1*sin(q1 + pi/4)*sin(q2)*sin(q4)*cos(q3) - L4*dq1*sin(q1 + pi/4)*cos(q2)*cos(q4) - L4*dq1*sin(q3)*sin(q4)*cos(q1 + pi/4) - L4*dq2*sin(q2)*cos(q1 + pi/4)*cos(q4) - L4*dq2*sin(q4)*cos(q1 + pi/4)*cos(q2)*cos(q3) - L4*dq3*sin(q1 + pi/4)*sin(q4)*cos(q3) + L4*dq3*sin(q2)*sin(q3)*sin(q4)*cos(q1 + pi/4) - L4*dq4*sin(q1 + pi/4)*sin(q3)*cos(q4) - L4*dq4*sin(q2)*cos(q1 + pi/4)*cos(q3)*cos(q4) - L4*dq4*sin(q4)*cos(q1 + pi/4)*cos(q2) - L5*dq1*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q5) + L5*dq1*sin(q1 + pi/4)*sin(q2)*cos(q3)*cos(q4)*cos(q5) + L5*dq1*sin(q1 + pi/4)*sin(q4)*cos(q2)*cos(q5) - L5*dq1*sin(q3)*cos(q1 + pi/4)*cos(q4)*cos(q5) - L5*dq1*sin(q5)*cos(q1 + pi/4)*cos(q3) + L5*dq2*sin(q2)*sin(q4)*cos(q1 + pi/4)*cos(q5) + L5*dq2*sin(q3)*sin(q5)*cos(q1 + pi/4)*cos(q2) - L5*dq2*cos(q1 + pi/4)*cos(q2)*cos(q3)*cos(q4)*cos(q5) + L5*dq3*sin(q1 + pi/4)*sin(q3)*sin(q5) - L5*dq3*sin(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) + L5*dq3*sin(q2)*sin(q3)*cos(q1 + pi/4)*cos(q4)*cos(q5) + L5*dq3*sin(q2)*sin(q5)*cos(q1 + pi/4)*cos(q3) + L5*dq4*sin(q1 + pi/4)*sin(q3)*sin(q4)*cos(q5) + L5*dq4*sin(q2)*sin(q4)*cos(q1 + pi/4)*cos(q3)*cos(q5) - L5*dq4*cos(q1 + pi/4)*cos(q2)*cos(q4)*cos(q5) + L5*dq5*sin(q1 + pi/4)*sin(q3)*sin(q5)*cos(q4) - L5*dq5*sin(q1 + pi/4)*cos(q3)*cos(q5) + L5*dq5*sin(q2)*sin(q3)*cos(q1 + pi/4)*cos(q5) + L5*dq5*sin(q2)*sin(q5)*cos(q1 + pi/4)*cos(q3)*cos(q4) + L5*dq5*sin(q4)*sin(q5)*cos(q1 + pi/4)*cos(q2), 
                       dq1*(-L2*sin(q2) - L3*cos(q2)*cos(q3) - L4*sin(q2)*cos(q4) - L4*sin(q4)*cos(q2)*cos(q3) + L5*sin(q2)*sin(q4)*cos(q5) + L5*sin(q3)*sin(q5)*cos(q2) - L5*cos(q2)*cos(q3)*cos(q4)*cos(q5))*cos(q1 + pi/4) + (-L2*dq2*cos(q2) + L3*dq2*sin(q2)*cos(q3) + L3*dq3*sin(q3)*cos(q2) + L4*dq2*sin(q2)*sin(q4)*cos(q3) - L4*dq2*cos(q2)*cos(q4) + L4*dq3*sin(q3)*sin(q4)*cos(q2) + L4*dq4*sin(q2)*sin(q4) - L4*dq4*cos(q2)*cos(q3)*cos(q4) - L5*dq2*sin(q2)*sin(q3)*sin(q5) + L5*dq2*sin(q2)*cos(q3)*cos(q4)*cos(q5) + L5*dq2*sin(q4)*cos(q2)*cos(q5) + L5*dq3*sin(q3)*cos(q2)*cos(q4)*cos(q5) + L5*dq3*sin(q5)*cos(q2)*cos(q3) + L5*dq4*sin(q2)*cos(q4)*cos(q5) + L5*dq4*sin(q4)*cos(q2)*cos(q3)*cos(q5) - L5*dq5*sin(q2)*sin(q4)*sin(q5) + L5*dq5*sin(q3)*cos(q2)*cos(q5) + L5*dq5*sin(q5)*cos(q2)*cos(q3)*cos(q4))*sin(q1 + pi/4), 
                       -L3*dq1*sin(q1 + pi/4)*cos(q3) + L3*dq1*sin(q2)*sin(q3)*cos(q1 + pi/4) + L3*dq2*sin(q1 + pi/4)*sin(q3)*cos(q2) + L3*dq3*sin(q1 + pi/4)*sin(q2)*cos(q3) - L3*dq3*sin(q3)*cos(q1 + pi/4) - L4*dq1*sin(q1 + pi/4)*sin(q4)*cos(q3) + L4*dq1*sin(q2)*sin(q3)*sin(q4)*cos(q1 + pi/4) + L4*dq2*sin(q1 + pi/4)*sin(q3)*sin(q4)*cos(q2) + L4*dq3*sin(q1 + pi/4)*sin(q2)*sin(q4)*cos(q3) - L4*dq3*sin(q3)*sin(q4)*cos(q1 + pi/4) + L4*dq4*sin(q1 + pi/4)*sin(q2)*sin(q3)*cos(q4) + L4*dq4*cos(q1 + pi/4)*cos(q3)*cos(q4) + L5*dq1*sin(q1 + pi/4)*sin(q3)*sin(q5) - L5*dq1*sin(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) + L5*dq1*sin(q2)*sin(q3)*cos(q1 + pi/4)*cos(q4)*cos(q5) + L5*dq1*sin(q2)*sin(q5)*cos(q1 + pi/4)*cos(q3) + L5*dq2*sin(q1 + pi/4)*sin(q3)*cos(q2)*cos(q4)*cos(q5) + L5*dq2*sin(q1 + pi/4)*sin(q5)*cos(q2)*cos(q3) - L5*dq3*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q5) + L5*dq3*sin(q1 + pi/4)*sin(q2)*cos(q3)*cos(q4)*cos(q5) - L5*dq3*sin(q3)*cos(q1 + pi/4)*cos(q4)*cos(q5) - L5*dq3*sin(q5)*cos(q1 + pi/4)*cos(q3) - L5*dq4*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q4)*cos(q5) - L5*dq4*sin(q4)*cos(q1 + pi/4)*cos(q3)*cos(q5) - L5*dq5*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q5)*cos(q4) + L5*dq5*sin(q1 + pi/4)*sin(q2)*cos(q3)*cos(q5) - L5*dq5*sin(q3)*cos(q1 + pi/4)*cos(q5) - L5*dq5*sin(q5)*cos(q1 + pi/4)*cos(q3)*cos(q4), 
                       -L4*dq1*sin(q1 + pi/4)*sin(q3)*cos(q4) - L4*dq1*sin(q2)*cos(q1 + pi/4)*cos(q3)*cos(q4) - L4*dq1*sin(q4)*cos(q1 + pi/4)*cos(q2) + L4*dq2*sin(q1 + pi/4)*sin(q2)*sin(q4) - L4*dq2*sin(q1 + pi/4)*cos(q2)*cos(q3)*cos(q4) + L4*dq3*sin(q1 + pi/4)*sin(q2)*sin(q3)*cos(q4) + L4*dq3*cos(q1 + pi/4)*cos(q3)*cos(q4) + L4*dq4*sin(q1 + pi/4)*sin(q2)*sin(q4)*cos(q3) - L4*dq4*sin(q1 + pi/4)*cos(q2)*cos(q4) - L4*dq4*sin(q3)*sin(q4)*cos(q1 + pi/4) + L5*dq1*sin(q1 + pi/4)*sin(q3)*sin(q4)*cos(q5) + L5*dq1*sin(q2)*sin(q4)*cos(q1 + pi/4)*cos(q3)*cos(q5) - L5*dq1*cos(q1 + pi/4)*cos(q2)*cos(q4)*cos(q5) + L5*dq2*sin(q1 + pi/4)*sin(q2)*cos(q4)*cos(q5) + L5*dq2*sin(q1 + pi/4)*sin(q4)*cos(q2)*cos(q3)*cos(q5) - L5*dq3*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q4)*cos(q5) - L5*dq3*sin(q4)*cos(q1 + pi/4)*cos(q3)*cos(q5) + L5*dq4*sin(q1 + pi/4)*sin(q2)*cos(q3)*cos(q4)*cos(q5) + L5*dq4*sin(q1 + pi/4)*sin(q4)*cos(q2)*cos(q5) - L5*dq4*sin(q3)*cos(q1 + pi/4)*cos(q4)*cos(q5) - L5*dq5*sin(q1 + pi/4)*sin(q2)*sin(q4)*sin(q5)*cos(q3) + L5*dq5*sin(q1 + pi/4)*sin(q5)*cos(q2)*cos(q4) + L5*dq5*sin(q3)*sin(q4)*sin(q5)*cos(q1 + pi/4), 
                       L5*(dq1*sin(q1 + pi/4)*sin(q3)*sin(q5)*cos(q4) - dq1*sin(q1 + pi/4)*cos(q3)*cos(q5) + dq1*sin(q2)*sin(q3)*cos(q1 + pi/4)*cos(q5) + dq1*sin(q2)*sin(q5)*cos(q1 + pi/4)*cos(q3)*cos(q4) + dq1*sin(q4)*sin(q5)*cos(q1 + pi/4)*cos(q2) - dq2*sin(q1 + pi/4)*sin(q2)*sin(q4)*sin(q5) + dq2*sin(q1 + pi/4)*sin(q3)*cos(q2)*cos(q5) + dq2*sin(q1 + pi/4)*sin(q5)*cos(q2)*cos(q3)*cos(q4) - dq3*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q5)*cos(q4) + dq3*sin(q1 + pi/4)*sin(q2)*cos(q3)*cos(q5) - dq3*sin(q3)*cos(q1 + pi/4)*cos(q5) - dq3*sin(q5)*cos(q1 + pi/4)*cos(q3)*cos(q4) - dq4*sin(q1 + pi/4)*sin(q2)*sin(q4)*sin(q5)*cos(q3) + dq4*sin(q1 + pi/4)*sin(q5)*cos(q2)*cos(q4) + dq4*sin(q3)*sin(q4)*sin(q5)*cos(q1 + pi/4) - dq5*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q5) + dq5*sin(q1 + pi/4)*sin(q2)*cos(q3)*cos(q4)*cos(q5) + dq5*sin(q1 + pi/4)*sin(q4)*cos(q2)*cos(q5) - dq5*sin(q3)*cos(q1 + pi/4)*cos(q4)*cos(q5) - dq5*sin(q5)*cos(q1 + pi/4)*cos(q3)), 
                       0, 
                       0], 
                      [0, 
                       L2*dq2*sin(q2) + L3*dq2*cos(q2)*cos(q3) - L3*dq3*sin(q2)*sin(q3) + L4*(dq2*sin(q2)*cos(q4) + dq2*sin(q4)*cos(q2)*cos(q3) - dq3*sin(q2)*sin(q3)*sin(q4) + dq4*sin(q2)*cos(q3)*cos(q4) + dq4*sin(q4)*cos(q2)) + L5*(-dq2*sin(q3)*sin(q5)*cos(q2) - dq3*sin(q2)*sin(q5)*cos(q3) - dq5*(sin(q2)*cos(q3)*cos(q4) + sin(q4)*cos(q2))*sin(q5) - dq5*sin(q2)*sin(q3)*cos(q5) + (-dq2*sin(q2)*sin(q4) + dq2*cos(q2)*cos(q3)*cos(q4) - dq3*sin(q2)*sin(q3)*cos(q4) - dq4*sin(q2)*sin(q4)*cos(q3) + dq4*cos(q2)*cos(q4))*cos(q5)), 
                       -dq2*(L3*sin(q3) + L4*sin(q3)*sin(q4) + L5*(sin(q3)*cos(q4)*cos(q5) + sin(q5)*cos(q3)))*sin(q2) + (L3*dq3*cos(q3) + L4*dq3*sin(q4)*cos(q3) + L4*dq4*sin(q3)*cos(q4) + L5*(-dq3*sin(q3)*sin(q5) + dq3*cos(q3)*cos(q4)*cos(q5) - dq4*sin(q3)*sin(q4)*cos(q5) - dq5*sin(q3)*sin(q5)*cos(q4) + dq5*cos(q3)*cos(q5)))*cos(q2), 
                       L4*(dq2*sin(q2)*cos(q3)*cos(q4) + dq2*sin(q4)*cos(q2) + dq3*sin(q3)*cos(q2)*cos(q4) + dq4*sin(q2)*cos(q4) + dq4*sin(q4)*cos(q2)*cos(q3)) - L5*dq5*(sin(q2)*cos(q4) + sin(q4)*cos(q2)*cos(q3))*sin(q5) + L5*(-dq2*sin(q2)*sin(q4)*cos(q3) + dq2*cos(q2)*cos(q4) - dq3*sin(q3)*sin(q4)*cos(q2) - dq4*sin(q2)*sin(q4) + dq4*cos(q2)*cos(q3)*cos(q4))*cos(q5), -L5*(dq2*sin(q2)*sin(q3)*cos(q5) - dq3*cos(q2)*cos(q3)*cos(q5) + dq5*(sin(q2)*sin(q4) - cos(q2)*cos(q3)*cos(q4))*cos(q5) + dq5*sin(q3)*sin(q5)*cos(q2) + (dq2*sin(q2)*cos(q3)*cos(q4) + dq2*sin(q4)*cos(q2) + dq3*sin(q3)*cos(q2)*cos(q4) + dq4*sin(q2)*cos(q4) + dq4*sin(q4)*cos(q2)*cos(q3))*sin(q5)), 
                       0, 
                       0]])
        return z
    
    def djacobi_7(self, q, dq):
        L, h, H, L0, L1, L2, L3, L4, L5, L6 = self.L, self.h, self.H, self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6
        q1, q2, q3, q4, q5, q6, q7 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0], q[6, 0]
        dq1, dq2, dq3, dq4, dq5, dq6, dq7 = dq[0, 0], dq[1, 0], dq[2, 0], dq[3, 0], dq[4, 0], dq[5, 0], dq[6, 0]
        z = np.array([[-L1*dq1*cos(q1 + pi/4) - L2*dq1*cos(q1 + pi/4)*cos(q2) + L2*dq2*sin(q1 + pi/4)*sin(q2) + L3*dq1*sin(q1 + pi/4)*sin(q3) + L3*dq1*sin(q2)*cos(q1 + pi/4)*cos(q3) + L3*dq2*sin(q1 + pi/4)*cos(q2)*cos(q3) - L3*dq3*sin(q1 + pi/4)*sin(q2)*sin(q3) - L3*dq3*cos(q1 + pi/4)*cos(q3) + L4*dq1*sin(q1 + pi/4)*sin(q3)*sin(q4) + L4*dq1*sin(q2)*sin(q4)*cos(q1 + pi/4)*cos(q3) - L4*dq1*cos(q1 + pi/4)*cos(q2)*cos(q4) + L4*dq2*sin(q1 + pi/4)*sin(q2)*cos(q4) + L4*dq2*sin(q1 + pi/4)*sin(q4)*cos(q2)*cos(q3) - L4*dq3*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q4) - L4*dq3*sin(q4)*cos(q1 + pi/4)*cos(q3) + L4*dq4*sin(q1 + pi/4)*sin(q2)*cos(q3)*cos(q4) + L4*dq4*sin(q1 + pi/4)*sin(q4)*cos(q2) - L4*dq4*sin(q3)*cos(q1 + pi/4)*cos(q4) + L5*dq1*sin(q1 + pi/4)*sin(q3)*cos(q4)*cos(q5) + L5*dq1*sin(q1 + pi/4)*sin(q5)*cos(q3) - L5*dq1*sin(q2)*sin(q3)*sin(q5)*cos(q1 + pi/4) + L5*dq1*sin(q2)*cos(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) + L5*dq1*sin(q4)*cos(q1 + pi/4)*cos(q2)*cos(q5) - L5*dq2*sin(q1 + pi/4)*sin(q2)*sin(q4)*cos(q5) - L5*dq2*sin(q1 + pi/4)*sin(q3)*sin(q5)*cos(q2) + L5*dq2*sin(q1 + pi/4)*cos(q2)*cos(q3)*cos(q4)*cos(q5) - L5*dq3*sin(q1 + pi/4)*sin(q2)*sin(q3)*cos(q4)*cos(q5) - L5*dq3*sin(q1 + pi/4)*sin(q2)*sin(q5)*cos(q3) + L5*dq3*sin(q3)*sin(q5)*cos(q1 + pi/4) - L5*dq3*cos(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) - L5*dq4*sin(q1 + pi/4)*sin(q2)*sin(q4)*cos(q3)*cos(q5) + L5*dq4*sin(q1 + pi/4)*cos(q2)*cos(q4)*cos(q5) + L5*dq4*sin(q3)*sin(q4)*cos(q1 + pi/4)*cos(q5) - L5*dq5*sin(q1 + pi/4)*sin(q2)*sin(q3)*cos(q5) - L5*dq5*sin(q1 + pi/4)*sin(q2)*sin(q5)*cos(q3)*cos(q4) - L5*dq5*sin(q1 + pi/4)*sin(q4)*sin(q5)*cos(q2) + L5*dq5*sin(q3)*sin(q5)*cos(q1 + pi/4)*cos(q4) - L5*dq5*cos(q1 + pi/4)*cos(q3)*cos(q5), 
                       -dq1*(-L2*sin(q2) - L3*cos(q2)*cos(q3) - L4*sin(q2)*cos(q4) - L4*sin(q4)*cos(q2)*cos(q3) + L5*sin(q2)*sin(q4)*cos(q5) + L5*sin(q3)*sin(q5)*cos(q2) - L5*cos(q2)*cos(q3)*cos(q4)*cos(q5))*sin(q1 + pi/4) + (-L2*dq2*cos(q2) + L3*dq2*sin(q2)*cos(q3) + L3*dq3*sin(q3)*cos(q2) + L4*dq2*sin(q2)*sin(q4)*cos(q3) - L4*dq2*cos(q2)*cos(q4) + L4*dq3*sin(q3)*sin(q4)*cos(q2) + L4*dq4*sin(q2)*sin(q4) - L4*dq4*cos(q2)*cos(q3)*cos(q4) - L5*dq2*sin(q2)*sin(q3)*sin(q5) + L5*dq2*sin(q2)*cos(q3)*cos(q4)*cos(q5) + L5*dq2*sin(q4)*cos(q2)*cos(q5) + L5*dq3*sin(q3)*cos(q2)*cos(q4)*cos(q5) + L5*dq3*sin(q5)*cos(q2)*cos(q3) + L5*dq4*sin(q2)*cos(q4)*cos(q5) + L5*dq4*sin(q4)*cos(q2)*cos(q3)*cos(q5) - L5*dq5*sin(q2)*sin(q4)*sin(q5) + L5*dq5*sin(q3)*cos(q2)*cos(q5) + L5*dq5*sin(q5)*cos(q2)*cos(q3)*cos(q4))*cos(q1 + pi/4), 
                       -L3*dq1*sin(q1 + pi/4)*sin(q2)*sin(q3) - L3*dq1*cos(q1 + pi/4)*cos(q3) + L3*dq2*sin(q3)*cos(q1 + pi/4)*cos(q2) + L3*dq3*sin(q1 + pi/4)*sin(q3) + L3*dq3*sin(q2)*cos(q1 + pi/4)*cos(q3) - L4*dq1*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q4) - L4*dq1*sin(q4)*cos(q1 + pi/4)*cos(q3) + L4*dq2*sin(q3)*sin(q4)*cos(q1 + pi/4)*cos(q2) + L4*dq3*sin(q1 + pi/4)*sin(q3)*sin(q4) + L4*dq3*sin(q2)*sin(q4)*cos(q1 + pi/4)*cos(q3) - L4*dq4*sin(q1 + pi/4)*cos(q3)*cos(q4) + L4*dq4*sin(q2)*sin(q3)*cos(q1 + pi/4)*cos(q4) - L5*dq1*sin(q1 + pi/4)*sin(q2)*sin(q3)*cos(q4)*cos(q5) - L5*dq1*sin(q1 + pi/4)*sin(q2)*sin(q5)*cos(q3) + L5*dq1*sin(q3)*sin(q5)*cos(q1 + pi/4) - L5*dq1*cos(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) + L5*dq2*sin(q3)*cos(q1 + pi/4)*cos(q2)*cos(q4)*cos(q5) + L5*dq2*sin(q5)*cos(q1 + pi/4)*cos(q2)*cos(q3) + L5*dq3*sin(q1 + pi/4)*sin(q3)*cos(q4)*cos(q5) + L5*dq3*sin(q1 + pi/4)*sin(q5)*cos(q3) - L5*dq3*sin(q2)*sin(q3)*sin(q5)*cos(q1 + pi/4) + L5*dq3*sin(q2)*cos(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) + L5*dq4*sin(q1 + pi/4)*sin(q4)*cos(q3)*cos(q5) - L5*dq4*sin(q2)*sin(q3)*sin(q4)*cos(q1 + pi/4)*cos(q5) + L5*dq5*sin(q1 + pi/4)*sin(q3)*cos(q5) + L5*dq5*sin(q1 + pi/4)*sin(q5)*cos(q3)*cos(q4) - L5*dq5*sin(q2)*sin(q3)*sin(q5)*cos(q1 + pi/4)*cos(q4) + L5*dq5*sin(q2)*cos(q1 + pi/4)*cos(q3)*cos(q5), L4*dq1*sin(q1 + pi/4)*sin(q2)*cos(q3)*cos(q4) + L4*dq1*sin(q1 + pi/4)*sin(q4)*cos(q2) - L4*dq1*sin(q3)*cos(q1 + pi/4)*cos(q4) + L4*dq2*sin(q2)*sin(q4)*cos(q1 + pi/4) - L4*dq2*cos(q1 + pi/4)*cos(q2)*cos(q3)*cos(q4) - L4*dq3*sin(q1 + pi/4)*cos(q3)*cos(q4) + L4*dq3*sin(q2)*sin(q3)*cos(q1 + pi/4)*cos(q4) + L4*dq4*sin(q1 + pi/4)*sin(q3)*sin(q4) + L4*dq4*sin(q2)*sin(q4)*cos(q1 + pi/4)*cos(q3) - L4*dq4*cos(q1 + pi/4)*cos(q2)*cos(q4) - L5*dq1*sin(q1 + pi/4)*sin(q2)*sin(q4)*cos(q3)*cos(q5) + L5*dq1*sin(q1 + pi/4)*cos(q2)*cos(q4)*cos(q5) + L5*dq1*sin(q3)*sin(q4)*cos(q1 + pi/4)*cos(q5) + L5*dq2*sin(q2)*cos(q1 + pi/4)*cos(q4)*cos(q5) + L5*dq2*sin(q4)*cos(q1 + pi/4)*cos(q2)*cos(q3)*cos(q5) + L5*dq3*sin(q1 + pi/4)*sin(q4)*cos(q3)*cos(q5) - L5*dq3*sin(q2)*sin(q3)*sin(q4)*cos(q1 + pi/4)*cos(q5) + L5*dq4*sin(q1 + pi/4)*sin(q3)*cos(q4)*cos(q5) + L5*dq4*sin(q2)*cos(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) + L5*dq4*sin(q4)*cos(q1 + pi/4)*cos(q2)*cos(q5) - L5*dq5*sin(q1 + pi/4)*sin(q3)*sin(q4)*sin(q5) - L5*dq5*sin(q2)*sin(q4)*sin(q5)*cos(q1 + pi/4)*cos(q3) + L5*dq5*sin(q5)*cos(q1 + pi/4)*cos(q2)*cos(q4), 
                       L5*(-dq1*sin(q1 + pi/4)*sin(q2)*sin(q3)*cos(q5) - dq1*sin(q1 + pi/4)*sin(q2)*sin(q5)*cos(q3)*cos(q4) - dq1*sin(q1 + pi/4)*sin(q4)*sin(q5)*cos(q2) + dq1*sin(q3)*sin(q5)*cos(q1 + pi/4)*cos(q4) - dq1*cos(q1 + pi/4)*cos(q3)*cos(q5) - dq2*sin(q2)*sin(q4)*sin(q5)*cos(q1 + pi/4) + dq2*sin(q3)*cos(q1 + pi/4)*cos(q2)*cos(q5) + dq2*sin(q5)*cos(q1 + pi/4)*cos(q2)*cos(q3)*cos(q4) + dq3*sin(q1 + pi/4)*sin(q3)*cos(q5) + dq3*sin(q1 + pi/4)*sin(q5)*cos(q3)*cos(q4) - dq3*sin(q2)*sin(q3)*sin(q5)*cos(q1 + pi/4)*cos(q4) + dq3*sin(q2)*cos(q1 + pi/4)*cos(q3)*cos(q5) - dq4*sin(q1 + pi/4)*sin(q3)*sin(q4)*sin(q5) - dq4*sin(q2)*sin(q4)*sin(q5)*cos(q1 + pi/4)*cos(q3) + dq4*sin(q5)*cos(q1 + pi/4)*cos(q2)*cos(q4) + dq5*sin(q1 + pi/4)*sin(q3)*cos(q4)*cos(q5) + dq5*sin(q1 + pi/4)*sin(q5)*cos(q3) - dq5*sin(q2)*sin(q3)*sin(q5)*cos(q1 + pi/4) + dq5*sin(q2)*cos(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) + dq5*sin(q4)*cos(q1 + pi/4)*cos(q2)*cos(q5)), 
                       0, 
                       0], 
                      [-L1*dq1*sin(q1 + pi/4) - L2*dq1*sin(q1 + pi/4)*cos(q2) - L2*dq2*sin(q2)*cos(q1 + pi/4) + L3*dq1*sin(q1 + pi/4)*sin(q2)*cos(q3) - L3*dq1*sin(q3)*cos(q1 + pi/4) - L3*dq2*cos(q1 + pi/4)*cos(q2)*cos(q3) - L3*dq3*sin(q1 + pi/4)*cos(q3) + L3*dq3*sin(q2)*sin(q3)*cos(q1 + pi/4) + L4*dq1*sin(q1 + pi/4)*sin(q2)*sin(q4)*cos(q3) - L4*dq1*sin(q1 + pi/4)*cos(q2)*cos(q4) - L4*dq1*sin(q3)*sin(q4)*cos(q1 + pi/4) - L4*dq2*sin(q2)*cos(q1 + pi/4)*cos(q4) - L4*dq2*sin(q4)*cos(q1 + pi/4)*cos(q2)*cos(q3) - L4*dq3*sin(q1 + pi/4)*sin(q4)*cos(q3) + L4*dq3*sin(q2)*sin(q3)*sin(q4)*cos(q1 + pi/4) - L4*dq4*sin(q1 + pi/4)*sin(q3)*cos(q4) - L4*dq4*sin(q2)*cos(q1 + pi/4)*cos(q3)*cos(q4) - L4*dq4*sin(q4)*cos(q1 + pi/4)*cos(q2) - L5*dq1*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q5) + L5*dq1*sin(q1 + pi/4)*sin(q2)*cos(q3)*cos(q4)*cos(q5) + L5*dq1*sin(q1 + pi/4)*sin(q4)*cos(q2)*cos(q5) - L5*dq1*sin(q3)*cos(q1 + pi/4)*cos(q4)*cos(q5) - L5*dq1*sin(q5)*cos(q1 + pi/4)*cos(q3) + L5*dq2*sin(q2)*sin(q4)*cos(q1 + pi/4)*cos(q5) + L5*dq2*sin(q3)*sin(q5)*cos(q1 + pi/4)*cos(q2) - L5*dq2*cos(q1 + pi/4)*cos(q2)*cos(q3)*cos(q4)*cos(q5) + L5*dq3*sin(q1 + pi/4)*sin(q3)*sin(q5) - L5*dq3*sin(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) + L5*dq3*sin(q2)*sin(q3)*cos(q1 + pi/4)*cos(q4)*cos(q5) + L5*dq3*sin(q2)*sin(q5)*cos(q1 + pi/4)*cos(q3) + L5*dq4*sin(q1 + pi/4)*sin(q3)*sin(q4)*cos(q5) + L5*dq4*sin(q2)*sin(q4)*cos(q1 + pi/4)*cos(q3)*cos(q5) - L5*dq4*cos(q1 + pi/4)*cos(q2)*cos(q4)*cos(q5) + L5*dq5*sin(q1 + pi/4)*sin(q3)*sin(q5)*cos(q4) - L5*dq5*sin(q1 + pi/4)*cos(q3)*cos(q5) + L5*dq5*sin(q2)*sin(q3)*cos(q1 + pi/4)*cos(q5) + L5*dq5*sin(q2)*sin(q5)*cos(q1 + pi/4)*cos(q3)*cos(q4) + L5*dq5*sin(q4)*sin(q5)*cos(q1 + pi/4)*cos(q2), 
                       dq1*(-L2*sin(q2) - L3*cos(q2)*cos(q3) - L4*sin(q2)*cos(q4) - L4*sin(q4)*cos(q2)*cos(q3) + L5*sin(q2)*sin(q4)*cos(q5) + L5*sin(q3)*sin(q5)*cos(q2) - L5*cos(q2)*cos(q3)*cos(q4)*cos(q5))*cos(q1 + pi/4) + (-L2*dq2*cos(q2) + L3*dq2*sin(q2)*cos(q3) + L3*dq3*sin(q3)*cos(q2) + L4*dq2*sin(q2)*sin(q4)*cos(q3) - L4*dq2*cos(q2)*cos(q4) + L4*dq3*sin(q3)*sin(q4)*cos(q2) + L4*dq4*sin(q2)*sin(q4) - L4*dq4*cos(q2)*cos(q3)*cos(q4) - L5*dq2*sin(q2)*sin(q3)*sin(q5) + L5*dq2*sin(q2)*cos(q3)*cos(q4)*cos(q5) + L5*dq2*sin(q4)*cos(q2)*cos(q5) + L5*dq3*sin(q3)*cos(q2)*cos(q4)*cos(q5) + L5*dq3*sin(q5)*cos(q2)*cos(q3) + L5*dq4*sin(q2)*cos(q4)*cos(q5) + L5*dq4*sin(q4)*cos(q2)*cos(q3)*cos(q5) - L5*dq5*sin(q2)*sin(q4)*sin(q5) + L5*dq5*sin(q3)*cos(q2)*cos(q5) + L5*dq5*sin(q5)*cos(q2)*cos(q3)*cos(q4))*sin(q1 + pi/4), 
                       -L3*dq1*sin(q1 + pi/4)*cos(q3) + L3*dq1*sin(q2)*sin(q3)*cos(q1 + pi/4) + L3*dq2*sin(q1 + pi/4)*sin(q3)*cos(q2) + L3*dq3*sin(q1 + pi/4)*sin(q2)*cos(q3) - L3*dq3*sin(q3)*cos(q1 + pi/4) - L4*dq1*sin(q1 + pi/4)*sin(q4)*cos(q3) + L4*dq1*sin(q2)*sin(q3)*sin(q4)*cos(q1 + pi/4) + L4*dq2*sin(q1 + pi/4)*sin(q3)*sin(q4)*cos(q2) + L4*dq3*sin(q1 + pi/4)*sin(q2)*sin(q4)*cos(q3) - L4*dq3*sin(q3)*sin(q4)*cos(q1 + pi/4) + L4*dq4*sin(q1 + pi/4)*sin(q2)*sin(q3)*cos(q4) + L4*dq4*cos(q1 + pi/4)*cos(q3)*cos(q4) + L5*dq1*sin(q1 + pi/4)*sin(q3)*sin(q5) - L5*dq1*sin(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) + L5*dq1*sin(q2)*sin(q3)*cos(q1 + pi/4)*cos(q4)*cos(q5) + L5*dq1*sin(q2)*sin(q5)*cos(q1 + pi/4)*cos(q3) + L5*dq2*sin(q1 + pi/4)*sin(q3)*cos(q2)*cos(q4)*cos(q5) + L5*dq2*sin(q1 + pi/4)*sin(q5)*cos(q2)*cos(q3) - L5*dq3*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q5) + L5*dq3*sin(q1 + pi/4)*sin(q2)*cos(q3)*cos(q4)*cos(q5) - L5*dq3*sin(q3)*cos(q1 + pi/4)*cos(q4)*cos(q5) - L5*dq3*sin(q5)*cos(q1 + pi/4)*cos(q3) - L5*dq4*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q4)*cos(q5) - L5*dq4*sin(q4)*cos(q1 + pi/4)*cos(q3)*cos(q5) - L5*dq5*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q5)*cos(q4) + L5*dq5*sin(q1 + pi/4)*sin(q2)*cos(q3)*cos(q5) - L5*dq5*sin(q3)*cos(q1 + pi/4)*cos(q5) - L5*dq5*sin(q5)*cos(q1 + pi/4)*cos(q3)*cos(q4), -L4*dq1*sin(q1 + pi/4)*sin(q3)*cos(q4) - L4*dq1*sin(q2)*cos(q1 + pi/4)*cos(q3)*cos(q4) - L4*dq1*sin(q4)*cos(q1 + pi/4)*cos(q2) + L4*dq2*sin(q1 + pi/4)*sin(q2)*sin(q4) - L4*dq2*sin(q1 + pi/4)*cos(q2)*cos(q3)*cos(q4) + L4*dq3*sin(q1 + pi/4)*sin(q2)*sin(q3)*cos(q4) + L4*dq3*cos(q1 + pi/4)*cos(q3)*cos(q4) + L4*dq4*sin(q1 + pi/4)*sin(q2)*sin(q4)*cos(q3) - L4*dq4*sin(q1 + pi/4)*cos(q2)*cos(q4) - L4*dq4*sin(q3)*sin(q4)*cos(q1 + pi/4) + L5*dq1*sin(q1 + pi/4)*sin(q3)*sin(q4)*cos(q5) + L5*dq1*sin(q2)*sin(q4)*cos(q1 + pi/4)*cos(q3)*cos(q5) - L5*dq1*cos(q1 + pi/4)*cos(q2)*cos(q4)*cos(q5) + L5*dq2*sin(q1 + pi/4)*sin(q2)*cos(q4)*cos(q5) + L5*dq2*sin(q1 + pi/4)*sin(q4)*cos(q2)*cos(q3)*cos(q5) - L5*dq3*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q4)*cos(q5) - L5*dq3*sin(q4)*cos(q1 + pi/4)*cos(q3)*cos(q5) + L5*dq4*sin(q1 + pi/4)*sin(q2)*cos(q3)*cos(q4)*cos(q5) + L5*dq4*sin(q1 + pi/4)*sin(q4)*cos(q2)*cos(q5) - L5*dq4*sin(q3)*cos(q1 + pi/4)*cos(q4)*cos(q5) - L5*dq5*sin(q1 + pi/4)*sin(q2)*sin(q4)*sin(q5)*cos(q3) + L5*dq5*sin(q1 + pi/4)*sin(q5)*cos(q2)*cos(q4) + L5*dq5*sin(q3)*sin(q4)*sin(q5)*cos(q1 + pi/4), 
                       L5*(dq1*sin(q1 + pi/4)*sin(q3)*sin(q5)*cos(q4) - dq1*sin(q1 + pi/4)*cos(q3)*cos(q5) + dq1*sin(q2)*sin(q3)*cos(q1 + pi/4)*cos(q5) + dq1*sin(q2)*sin(q5)*cos(q1 + pi/4)*cos(q3)*cos(q4) + dq1*sin(q4)*sin(q5)*cos(q1 + pi/4)*cos(q2) - dq2*sin(q1 + pi/4)*sin(q2)*sin(q4)*sin(q5) + dq2*sin(q1 + pi/4)*sin(q3)*cos(q2)*cos(q5) + dq2*sin(q1 + pi/4)*sin(q5)*cos(q2)*cos(q3)*cos(q4) - dq3*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q5)*cos(q4) + dq3*sin(q1 + pi/4)*sin(q2)*cos(q3)*cos(q5) - dq3*sin(q3)*cos(q1 + pi/4)*cos(q5) - dq3*sin(q5)*cos(q1 + pi/4)*cos(q3)*cos(q4) - dq4*sin(q1 + pi/4)*sin(q2)*sin(q4)*sin(q5)*cos(q3) + dq4*sin(q1 + pi/4)*sin(q5)*cos(q2)*cos(q4) + dq4*sin(q3)*sin(q4)*sin(q5)*cos(q1 + pi/4) - dq5*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q5) + dq5*sin(q1 + pi/4)*sin(q2)*cos(q3)*cos(q4)*cos(q5) + dq5*sin(q1 + pi/4)*sin(q4)*cos(q2)*cos(q5) - dq5*sin(q3)*cos(q1 + pi/4)*cos(q4)*cos(q5) - dq5*sin(q5)*cos(q1 + pi/4)*cos(q3)), 
                       0, 
                       0], 
                      [0, 
                       L2*dq2*sin(q2) + L3*dq2*cos(q2)*cos(q3) - L3*dq3*sin(q2)*sin(q3) + L4*(dq2*sin(q2)*cos(q4) + dq2*sin(q4)*cos(q2)*cos(q3) - dq3*sin(q2)*sin(q3)*sin(q4) + dq4*sin(q2)*cos(q3)*cos(q4) + dq4*sin(q4)*cos(q2)) + L5*(-dq2*sin(q3)*sin(q5)*cos(q2) - dq3*sin(q2)*sin(q5)*cos(q3) - dq5*(sin(q2)*cos(q3)*cos(q4) + sin(q4)*cos(q2))*sin(q5) - dq5*sin(q2)*sin(q3)*cos(q5) + (-dq2*sin(q2)*sin(q4) + dq2*cos(q2)*cos(q3)*cos(q4) - dq3*sin(q2)*sin(q3)*cos(q4) - dq4*sin(q2)*sin(q4)*cos(q3) + dq4*cos(q2)*cos(q4))*cos(q5)), 
                       -dq2*(L3*sin(q3) + L4*sin(q3)*sin(q4) + L5*(sin(q3)*cos(q4)*cos(q5) + sin(q5)*cos(q3)))*sin(q2) + (L3*dq3*cos(q3) + L4*dq3*sin(q4)*cos(q3) + L4*dq4*sin(q3)*cos(q4) + L5*(-dq3*sin(q3)*sin(q5) + dq3*cos(q3)*cos(q4)*cos(q5) - dq4*sin(q3)*sin(q4)*cos(q5) - dq5*sin(q3)*sin(q5)*cos(q4) + dq5*cos(q3)*cos(q5)))*cos(q2), 
                       L4*(dq2*sin(q2)*cos(q3)*cos(q4) + dq2*sin(q4)*cos(q2) + dq3*sin(q3)*cos(q2)*cos(q4) + dq4*sin(q2)*cos(q4) + dq4*sin(q4)*cos(q2)*cos(q3)) - L5*dq5*(sin(q2)*cos(q4) + sin(q4)*cos(q2)*cos(q3))*sin(q5) + L5*(-dq2*sin(q2)*sin(q4)*cos(q3) + dq2*cos(q2)*cos(q4) - dq3*sin(q3)*sin(q4)*cos(q2) - dq4*sin(q2)*sin(q4) + dq4*cos(q2)*cos(q3)*cos(q4))*cos(q5), 
                       -L5*(dq2*sin(q2)*sin(q3)*cos(q5) - dq3*cos(q2)*cos(q3)*cos(q5) + dq5*(sin(q2)*sin(q4) - cos(q2)*cos(q3)*cos(q4))*cos(q5) + dq5*sin(q3)*sin(q5)*cos(q2) + (dq2*sin(q2)*cos(q3)*cos(q4) + dq2*sin(q4)*cos(q2) + dq3*sin(q3)*cos(q2)*cos(q4) + dq4*sin(q2)*cos(q4) + dq4*sin(q4)*cos(q2)*cos(q3))*sin(q5)), 
                       0, 
                       0]])
        return z
    
    def djacobi_GL(self, q, dq):
        L, h, H, L0, L1, L2, L3, L4, L5, L6 = self.L, self.h, self.H, self.L0, self.L1, self.L2, self.L3, self.L4, self.L5, self.L6
        q1, q2, q3, q4, q5, q6, q7 = q[0, 0], q[1, 0], q[2, 0], q[3, 0], q[4, 0], q[5, 0], q[6, 0]
        dq1, dq2, dq3, dq4, dq5, dq6, dq7 = dq[0, 0], dq[1, 0], dq[2, 0], dq[3, 0], dq[4, 0], dq[5, 0], dq[6, 0]
        z = np.array([[-L1*dq1*cos(q1 + pi/4) - L2*dq1*cos(q1 + pi/4)*cos(q2) + L2*dq2*sin(q1 + pi/4)*sin(q2) + L3*dq1*sin(q1 + pi/4)*sin(q3) + L3*dq1*sin(q2)*cos(q1 + pi/4)*cos(q3) + L3*dq2*sin(q1 + pi/4)*cos(q2)*cos(q3) - L3*dq3*sin(q1 + pi/4)*sin(q2)*sin(q3) - L3*dq3*cos(q1 + pi/4)*cos(q3) + L4*dq1*sin(q1 + pi/4)*sin(q3)*sin(q4) + L4*dq1*sin(q2)*sin(q4)*cos(q1 + pi/4)*cos(q3) - L4*dq1*cos(q1 + pi/4)*cos(q2)*cos(q4) + L4*dq2*sin(q1 + pi/4)*sin(q2)*cos(q4) + L4*dq2*sin(q1 + pi/4)*sin(q4)*cos(q2)*cos(q3) - L4*dq3*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q4) - L4*dq3*sin(q4)*cos(q1 + pi/4)*cos(q3) + L4*dq4*sin(q1 + pi/4)*sin(q2)*cos(q3)*cos(q4) + L4*dq4*sin(q1 + pi/4)*sin(q4)*cos(q2) - L4*dq4*sin(q3)*cos(q1 + pi/4)*cos(q4) + L5*dq1*sin(q1 + pi/4)*sin(q3)*cos(q4)*cos(q5) + L5*dq1*sin(q1 + pi/4)*sin(q5)*cos(q3) - L5*dq1*sin(q2)*sin(q3)*sin(q5)*cos(q1 + pi/4) + L5*dq1*sin(q2)*cos(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) + L5*dq1*sin(q4)*cos(q1 + pi/4)*cos(q2)*cos(q5) - L5*dq2*sin(q1 + pi/4)*sin(q2)*sin(q4)*cos(q5) - L5*dq2*sin(q1 + pi/4)*sin(q3)*sin(q5)*cos(q2) + L5*dq2*sin(q1 + pi/4)*cos(q2)*cos(q3)*cos(q4)*cos(q5) - L5*dq3*sin(q1 + pi/4)*sin(q2)*sin(q3)*cos(q4)*cos(q5) - L5*dq3*sin(q1 + pi/4)*sin(q2)*sin(q5)*cos(q3) + L5*dq3*sin(q3)*sin(q5)*cos(q1 + pi/4) - L5*dq3*cos(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) - L5*dq4*sin(q1 + pi/4)*sin(q2)*sin(q4)*cos(q3)*cos(q5) + L5*dq4*sin(q1 + pi/4)*cos(q2)*cos(q4)*cos(q5) + L5*dq4*sin(q3)*sin(q4)*cos(q1 + pi/4)*cos(q5) - L5*dq5*sin(q1 + pi/4)*sin(q2)*sin(q3)*cos(q5) - L5*dq5*sin(q1 + pi/4)*sin(q2)*sin(q5)*cos(q3)*cos(q4) - L5*dq5*sin(q1 + pi/4)*sin(q4)*sin(q5)*cos(q2) + L5*dq5*sin(q3)*sin(q5)*cos(q1 + pi/4)*cos(q4) - L5*dq5*cos(q1 + pi/4)*cos(q3)*cos(q5) + L6*dq1*sin(q1 + pi/4)*sin(q3)*sin(q4)*cos(q6) + L6*dq1*sin(q1 + pi/4)*sin(q3)*sin(q6)*cos(q4)*cos(q5) + L6*dq1*sin(q1 + pi/4)*sin(q5)*sin(q6)*cos(q3) - L6*dq1*sin(q2)*sin(q3)*sin(q5)*sin(q6)*cos(q1 + pi/4) + L6*dq1*sin(q2)*sin(q4)*cos(q1 + pi/4)*cos(q3)*cos(q6) + L6*dq1*sin(q2)*sin(q6)*cos(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) + L6*dq1*sin(q4)*sin(q6)*cos(q1 + pi/4)*cos(q2)*cos(q5) - L6*dq1*cos(q1 + pi/4)*cos(q2)*cos(q4)*cos(q6) - L6*dq2*sin(q1 + pi/4)*sin(q2)*sin(q4)*sin(q6)*cos(q5) + L6*dq2*sin(q1 + pi/4)*sin(q2)*cos(q4)*cos(q6) - L6*dq2*sin(q1 + pi/4)*sin(q3)*sin(q5)*sin(q6)*cos(q2) + L6*dq2*sin(q1 + pi/4)*sin(q4)*cos(q2)*cos(q3)*cos(q6) + L6*dq2*sin(q1 + pi/4)*sin(q6)*cos(q2)*cos(q3)*cos(q4)*cos(q5) - L6*dq3*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q4)*cos(q6) - L6*dq3*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q6)*cos(q4)*cos(q5) - L6*dq3*sin(q1 + pi/4)*sin(q2)*sin(q5)*sin(q6)*cos(q3) + L6*dq3*sin(q3)*sin(q5)*sin(q6)*cos(q1 + pi/4) - L6*dq3*sin(q4)*cos(q1 + pi/4)*cos(q3)*cos(q6) - L6*dq3*sin(q6)*cos(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) - L6*dq4*sin(q1 + pi/4)*sin(q2)*sin(q4)*sin(q6)*cos(q3)*cos(q5) + L6*dq4*sin(q1 + pi/4)*sin(q2)*cos(q3)*cos(q4)*cos(q6) + L6*dq4*sin(q1 + pi/4)*sin(q4)*cos(q2)*cos(q6) + L6*dq4*sin(q1 + pi/4)*sin(q6)*cos(q2)*cos(q4)*cos(q5) + L6*dq4*sin(q3)*sin(q4)*sin(q6)*cos(q1 + pi/4)*cos(q5) - L6*dq4*sin(q3)*cos(q1 + pi/4)*cos(q4)*cos(q6) - L6*dq5*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q6)*cos(q5) - L6*dq5*sin(q1 + pi/4)*sin(q2)*sin(q5)*sin(q6)*cos(q3)*cos(q4) - L6*dq5*sin(q1 + pi/4)*sin(q4)*sin(q5)*sin(q6)*cos(q2) + L6*dq5*sin(q3)*sin(q5)*sin(q6)*cos(q1 + pi/4)*cos(q4) - L6*dq5*sin(q6)*cos(q1 + pi/4)*cos(q3)*cos(q5) - L6*dq6*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q5)*cos(q6) - L6*dq6*sin(q1 + pi/4)*sin(q2)*sin(q4)*sin(q6)*cos(q3) + L6*dq6*sin(q1 + pi/4)*sin(q2)*cos(q3)*cos(q4)*cos(q5)*cos(q6) + L6*dq6*sin(q1 + pi/4)*sin(q4)*cos(q2)*cos(q5)*cos(q6) + L6*dq6*sin(q1 + pi/4)*sin(q6)*cos(q2)*cos(q4) + L6*dq6*sin(q3)*sin(q4)*sin(q6)*cos(q1 + pi/4) - L6*dq6*sin(q3)*cos(q1 + pi/4)*cos(q4)*cos(q5)*cos(q6) - L6*dq6*sin(q5)*cos(q1 + pi/4)*cos(q3)*cos(q6), 
                       -dq1*(-L2*sin(q2) - L3*cos(q2)*cos(q3) - L4*sin(q2)*cos(q4) - L4*sin(q4)*cos(q2)*cos(q3) + L5*sin(q2)*sin(q4)*cos(q5) + L5*sin(q3)*sin(q5)*cos(q2) - L5*cos(q2)*cos(q3)*cos(q4)*cos(q5) + L6*sin(q2)*sin(q4)*sin(q6)*cos(q5) - L6*sin(q2)*cos(q4)*cos(q6) + L6*sin(q3)*sin(q5)*sin(q6)*cos(q2) - L6*sin(q4)*cos(q2)*cos(q3)*cos(q6) - L6*sin(q6)*cos(q2)*cos(q3)*cos(q4)*cos(q5))*sin(q1 + pi/4) + (-L2*dq2*cos(q2) + L3*dq2*sin(q2)*cos(q3) + L3*dq3*sin(q3)*cos(q2) + L4*dq2*sin(q2)*sin(q4)*cos(q3) - L4*dq2*cos(q2)*cos(q4) + L4*dq3*sin(q3)*sin(q4)*cos(q2) + L4*dq4*sin(q2)*sin(q4) - L4*dq4*cos(q2)*cos(q3)*cos(q4) - L5*dq2*sin(q2)*sin(q3)*sin(q5) + L5*dq2*sin(q2)*cos(q3)*cos(q4)*cos(q5) + L5*dq2*sin(q4)*cos(q2)*cos(q5) + L5*dq3*sin(q3)*cos(q2)*cos(q4)*cos(q5) + L5*dq3*sin(q5)*cos(q2)*cos(q3) + L5*dq4*sin(q2)*cos(q4)*cos(q5) + L5*dq4*sin(q4)*cos(q2)*cos(q3)*cos(q5) - L5*dq5*sin(q2)*sin(q4)*sin(q5) + L5*dq5*sin(q3)*cos(q2)*cos(q5) + L5*dq5*sin(q5)*cos(q2)*cos(q3)*cos(q4) - L6*dq2*sin(q2)*sin(q3)*sin(q5)*sin(q6) + L6*dq2*sin(q2)*sin(q4)*cos(q3)*cos(q6) + L6*dq2*sin(q2)*sin(q6)*cos(q3)*cos(q4)*cos(q5) + L6*dq2*sin(q4)*sin(q6)*cos(q2)*cos(q5) - L6*dq2*cos(q2)*cos(q4)*cos(q6) + L6*dq3*sin(q3)*sin(q4)*cos(q2)*cos(q6) + L6*dq3*sin(q3)*sin(q6)*cos(q2)*cos(q4)*cos(q5) + L6*dq3*sin(q5)*sin(q6)*cos(q2)*cos(q3) + L6*dq4*sin(q2)*sin(q4)*cos(q6) + L6*dq4*sin(q2)*sin(q6)*cos(q4)*cos(q5) + L6*dq4*sin(q4)*sin(q6)*cos(q2)*cos(q3)*cos(q5) - L6*dq4*cos(q2)*cos(q3)*cos(q4)*cos(q6) - L6*dq5*sin(q2)*sin(q4)*sin(q5)*sin(q6) + L6*dq5*sin(q3)*sin(q6)*cos(q2)*cos(q5) + L6*dq5*sin(q5)*sin(q6)*cos(q2)*cos(q3)*cos(q4) + L6*dq6*sin(q2)*sin(q4)*cos(q5)*cos(q6) + L6*dq6*sin(q2)*sin(q6)*cos(q4) + L6*dq6*sin(q3)*sin(q5)*cos(q2)*cos(q6) + L6*dq6*sin(q4)*sin(q6)*cos(q2)*cos(q3) - L6*dq6*cos(q2)*cos(q3)*cos(q4)*cos(q5)*cos(q6))*cos(q1 + pi/4), 
                       -L3*dq1*sin(q1 + pi/4)*sin(q2)*sin(q3) - L3*dq1*cos(q1 + pi/4)*cos(q3) + L3*dq2*sin(q3)*cos(q1 + pi/4)*cos(q2) + L3*dq3*sin(q1 + pi/4)*sin(q3) + L3*dq3*sin(q2)*cos(q1 + pi/4)*cos(q3) - L4*dq1*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q4) - L4*dq1*sin(q4)*cos(q1 + pi/4)*cos(q3) + L4*dq2*sin(q3)*sin(q4)*cos(q1 + pi/4)*cos(q2) + L4*dq3*sin(q1 + pi/4)*sin(q3)*sin(q4) + L4*dq3*sin(q2)*sin(q4)*cos(q1 + pi/4)*cos(q3) - L4*dq4*sin(q1 + pi/4)*cos(q3)*cos(q4) + L4*dq4*sin(q2)*sin(q3)*cos(q1 + pi/4)*cos(q4) - L5*dq1*sin(q1 + pi/4)*sin(q2)*sin(q3)*cos(q4)*cos(q5) - L5*dq1*sin(q1 + pi/4)*sin(q2)*sin(q5)*cos(q3) + L5*dq1*sin(q3)*sin(q5)*cos(q1 + pi/4) - L5*dq1*cos(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) + L5*dq2*sin(q3)*cos(q1 + pi/4)*cos(q2)*cos(q4)*cos(q5) + L5*dq2*sin(q5)*cos(q1 + pi/4)*cos(q2)*cos(q3) + L5*dq3*sin(q1 + pi/4)*sin(q3)*cos(q4)*cos(q5) + L5*dq3*sin(q1 + pi/4)*sin(q5)*cos(q3) - L5*dq3*sin(q2)*sin(q3)*sin(q5)*cos(q1 + pi/4) + L5*dq3*sin(q2)*cos(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) + L5*dq4*sin(q1 + pi/4)*sin(q4)*cos(q3)*cos(q5) - L5*dq4*sin(q2)*sin(q3)*sin(q4)*cos(q1 + pi/4)*cos(q5) + L5*dq5*sin(q1 + pi/4)*sin(q3)*cos(q5) + L5*dq5*sin(q1 + pi/4)*sin(q5)*cos(q3)*cos(q4) - L5*dq5*sin(q2)*sin(q3)*sin(q5)*cos(q1 + pi/4)*cos(q4) + L5*dq5*sin(q2)*cos(q1 + pi/4)*cos(q3)*cos(q5) - L6*dq1*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q4)*cos(q6) - L6*dq1*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q6)*cos(q4)*cos(q5) - L6*dq1*sin(q1 + pi/4)*sin(q2)*sin(q5)*sin(q6)*cos(q3) + L6*dq1*sin(q3)*sin(q5)*sin(q6)*cos(q1 + pi/4) - L6*dq1*sin(q4)*cos(q1 + pi/4)*cos(q3)*cos(q6) - L6*dq1*sin(q6)*cos(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) + L6*dq2*sin(q3)*sin(q4)*cos(q1 + pi/4)*cos(q2)*cos(q6) + L6*dq2*sin(q3)*sin(q6)*cos(q1 + pi/4)*cos(q2)*cos(q4)*cos(q5) + L6*dq2*sin(q5)*sin(q6)*cos(q1 + pi/4)*cos(q2)*cos(q3) + L6*dq3*sin(q1 + pi/4)*sin(q3)*sin(q4)*cos(q6) + L6*dq3*sin(q1 + pi/4)*sin(q3)*sin(q6)*cos(q4)*cos(q5) + L6*dq3*sin(q1 + pi/4)*sin(q5)*sin(q6)*cos(q3) - L6*dq3*sin(q2)*sin(q3)*sin(q5)*sin(q6)*cos(q1 + pi/4) + L6*dq3*sin(q2)*sin(q4)*cos(q1 + pi/4)*cos(q3)*cos(q6) + L6*dq3*sin(q2)*sin(q6)*cos(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) + L6*dq4*sin(q1 + pi/4)*sin(q4)*sin(q6)*cos(q3)*cos(q5) - L6*dq4*sin(q1 + pi/4)*cos(q3)*cos(q4)*cos(q6) - L6*dq4*sin(q2)*sin(q3)*sin(q4)*sin(q6)*cos(q1 + pi/4)*cos(q5) + L6*dq4*sin(q2)*sin(q3)*cos(q1 + pi/4)*cos(q4)*cos(q6) + L6*dq5*sin(q1 + pi/4)*sin(q3)*sin(q6)*cos(q5) + L6*dq5*sin(q1 + pi/4)*sin(q5)*sin(q6)*cos(q3)*cos(q4) - L6*dq5*sin(q2)*sin(q3)*sin(q5)*sin(q6)*cos(q1 + pi/4)*cos(q4) + L6*dq5*sin(q2)*sin(q6)*cos(q1 + pi/4)*cos(q3)*cos(q5) + L6*dq6*sin(q1 + pi/4)*sin(q3)*sin(q5)*cos(q6) + L6*dq6*sin(q1 + pi/4)*sin(q4)*sin(q6)*cos(q3) - L6*dq6*sin(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5)*cos(q6) - L6*dq6*sin(q2)*sin(q3)*sin(q4)*sin(q6)*cos(q1 + pi/4) + L6*dq6*sin(q2)*sin(q3)*cos(q1 + pi/4)*cos(q4)*cos(q5)*cos(q6) + L6*dq6*sin(q2)*sin(q5)*cos(q1 + pi/4)*cos(q3)*cos(q6), 
                       L4*dq1*sin(q1 + pi/4)*sin(q2)*cos(q3)*cos(q4) + L4*dq1*sin(q1 + pi/4)*sin(q4)*cos(q2) - L4*dq1*sin(q3)*cos(q1 + pi/4)*cos(q4) + L4*dq2*sin(q2)*sin(q4)*cos(q1 + pi/4) - L4*dq2*cos(q1 + pi/4)*cos(q2)*cos(q3)*cos(q4) - L4*dq3*sin(q1 + pi/4)*cos(q3)*cos(q4) + L4*dq3*sin(q2)*sin(q3)*cos(q1 + pi/4)*cos(q4) + L4*dq4*sin(q1 + pi/4)*sin(q3)*sin(q4) + L4*dq4*sin(q2)*sin(q4)*cos(q1 + pi/4)*cos(q3) - L4*dq4*cos(q1 + pi/4)*cos(q2)*cos(q4) - L5*dq1*sin(q1 + pi/4)*sin(q2)*sin(q4)*cos(q3)*cos(q5) + L5*dq1*sin(q1 + pi/4)*cos(q2)*cos(q4)*cos(q5) + L5*dq1*sin(q3)*sin(q4)*cos(q1 + pi/4)*cos(q5) + L5*dq2*sin(q2)*cos(q1 + pi/4)*cos(q4)*cos(q5) + L5*dq2*sin(q4)*cos(q1 + pi/4)*cos(q2)*cos(q3)*cos(q5) + L5*dq3*sin(q1 + pi/4)*sin(q4)*cos(q3)*cos(q5) - L5*dq3*sin(q2)*sin(q3)*sin(q4)*cos(q1 + pi/4)*cos(q5) + L5*dq4*sin(q1 + pi/4)*sin(q3)*cos(q4)*cos(q5) + L5*dq4*sin(q2)*cos(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) + L5*dq4*sin(q4)*cos(q1 + pi/4)*cos(q2)*cos(q5) - L5*dq5*sin(q1 + pi/4)*sin(q3)*sin(q4)*sin(q5) - L5*dq5*sin(q2)*sin(q4)*sin(q5)*cos(q1 + pi/4)*cos(q3) + L5*dq5*sin(q5)*cos(q1 + pi/4)*cos(q2)*cos(q4) - L6*dq1*sin(q1 + pi/4)*sin(q2)*sin(q4)*sin(q6)*cos(q3)*cos(q5) + L6*dq1*sin(q1 + pi/4)*sin(q2)*cos(q3)*cos(q4)*cos(q6) + L6*dq1*sin(q1 + pi/4)*sin(q4)*cos(q2)*cos(q6) + L6*dq1*sin(q1 + pi/4)*sin(q6)*cos(q2)*cos(q4)*cos(q5) + L6*dq1*sin(q3)*sin(q4)*sin(q6)*cos(q1 + pi/4)*cos(q5) - L6*dq1*sin(q3)*cos(q1 + pi/4)*cos(q4)*cos(q6) + L6*dq2*sin(q2)*sin(q4)*cos(q1 + pi/4)*cos(q6) + L6*dq2*sin(q2)*sin(q6)*cos(q1 + pi/4)*cos(q4)*cos(q5) + L6*dq2*sin(q4)*sin(q6)*cos(q1 + pi/4)*cos(q2)*cos(q3)*cos(q5) - L6*dq2*cos(q1 + pi/4)*cos(q2)*cos(q3)*cos(q4)*cos(q6) + L6*dq3*sin(q1 + pi/4)*sin(q4)*sin(q6)*cos(q3)*cos(q5) - L6*dq3*sin(q1 + pi/4)*cos(q3)*cos(q4)*cos(q6) - L6*dq3*sin(q2)*sin(q3)*sin(q4)*sin(q6)*cos(q1 + pi/4)*cos(q5) + L6*dq3*sin(q2)*sin(q3)*cos(q1 + pi/4)*cos(q4)*cos(q6) + L6*dq4*sin(q1 + pi/4)*sin(q3)*sin(q4)*cos(q6) + L6*dq4*sin(q1 + pi/4)*sin(q3)*sin(q6)*cos(q4)*cos(q5) + L6*dq4*sin(q2)*sin(q4)*cos(q1 + pi/4)*cos(q3)*cos(q6) + L6*dq4*sin(q2)*sin(q6)*cos(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) + L6*dq4*sin(q4)*sin(q6)*cos(q1 + pi/4)*cos(q2)*cos(q5) - L6*dq4*cos(q1 + pi/4)*cos(q2)*cos(q4)*cos(q6) - L6*dq5*sin(q1 + pi/4)*sin(q3)*sin(q4)*sin(q5)*sin(q6) - L6*dq5*sin(q2)*sin(q4)*sin(q5)*sin(q6)*cos(q1 + pi/4)*cos(q3) + L6*dq5*sin(q5)*sin(q6)*cos(q1 + pi/4)*cos(q2)*cos(q4) + L6*dq6*sin(q1 + pi/4)*sin(q3)*sin(q4)*cos(q5)*cos(q6) + L6*dq6*sin(q1 + pi/4)*sin(q3)*sin(q6)*cos(q4) + L6*dq6*sin(q2)*sin(q4)*cos(q1 + pi/4)*cos(q3)*cos(q5)*cos(q6) + L6*dq6*sin(q2)*sin(q6)*cos(q1 + pi/4)*cos(q3)*cos(q4) + L6*dq6*sin(q4)*sin(q6)*cos(q1 + pi/4)*cos(q2) - L6*dq6*cos(q1 + pi/4)*cos(q2)*cos(q4)*cos(q5)*cos(q6), 
                       -L5*dq1*sin(q1 + pi/4)*sin(q2)*sin(q3)*cos(q5) - L5*dq1*sin(q1 + pi/4)*sin(q2)*sin(q5)*cos(q3)*cos(q4) - L5*dq1*sin(q1 + pi/4)*sin(q4)*sin(q5)*cos(q2) + L5*dq1*sin(q3)*sin(q5)*cos(q1 + pi/4)*cos(q4) - L5*dq1*cos(q1 + pi/4)*cos(q3)*cos(q5) - L5*dq2*sin(q2)*sin(q4)*sin(q5)*cos(q1 + pi/4) + L5*dq2*sin(q3)*cos(q1 + pi/4)*cos(q2)*cos(q5) + L5*dq2*sin(q5)*cos(q1 + pi/4)*cos(q2)*cos(q3)*cos(q4) + L5*dq3*sin(q1 + pi/4)*sin(q3)*cos(q5) + L5*dq3*sin(q1 + pi/4)*sin(q5)*cos(q3)*cos(q4) - L5*dq3*sin(q2)*sin(q3)*sin(q5)*cos(q1 + pi/4)*cos(q4) + L5*dq3*sin(q2)*cos(q1 + pi/4)*cos(q3)*cos(q5) - L5*dq4*sin(q1 + pi/4)*sin(q3)*sin(q4)*sin(q5) - L5*dq4*sin(q2)*sin(q4)*sin(q5)*cos(q1 + pi/4)*cos(q3) + L5*dq4*sin(q5)*cos(q1 + pi/4)*cos(q2)*cos(q4) + L5*dq5*sin(q1 + pi/4)*sin(q3)*cos(q4)*cos(q5) + L5*dq5*sin(q1 + pi/4)*sin(q5)*cos(q3) - L5*dq5*sin(q2)*sin(q3)*sin(q5)*cos(q1 + pi/4) + L5*dq5*sin(q2)*cos(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) + L5*dq5*sin(q4)*cos(q1 + pi/4)*cos(q2)*cos(q5) - L6*dq1*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q6)*cos(q5) - L6*dq1*sin(q1 + pi/4)*sin(q2)*sin(q5)*sin(q6)*cos(q3)*cos(q4) - L6*dq1*sin(q1 + pi/4)*sin(q4)*sin(q5)*sin(q6)*cos(q2) + L6*dq1*sin(q3)*sin(q5)*sin(q6)*cos(q1 + pi/4)*cos(q4) - L6*dq1*sin(q6)*cos(q1 + pi/4)*cos(q3)*cos(q5) - L6*dq2*sin(q2)*sin(q4)*sin(q5)*sin(q6)*cos(q1 + pi/4) + L6*dq2*sin(q3)*sin(q6)*cos(q1 + pi/4)*cos(q2)*cos(q5) + L6*dq2*sin(q5)*sin(q6)*cos(q1 + pi/4)*cos(q2)*cos(q3)*cos(q4) + L6*dq3*sin(q1 + pi/4)*sin(q3)*sin(q6)*cos(q5) + L6*dq3*sin(q1 + pi/4)*sin(q5)*sin(q6)*cos(q3)*cos(q4) - L6*dq3*sin(q2)*sin(q3)*sin(q5)*sin(q6)*cos(q1 + pi/4)*cos(q4) + L6*dq3*sin(q2)*sin(q6)*cos(q1 + pi/4)*cos(q3)*cos(q5) - L6*dq4*sin(q1 + pi/4)*sin(q3)*sin(q4)*sin(q5)*sin(q6) - L6*dq4*sin(q2)*sin(q4)*sin(q5)*sin(q6)*cos(q1 + pi/4)*cos(q3) + L6*dq4*sin(q5)*sin(q6)*cos(q1 + pi/4)*cos(q2)*cos(q4) + L6*dq5*sin(q1 + pi/4)*sin(q3)*sin(q6)*cos(q4)*cos(q5) + L6*dq5*sin(q1 + pi/4)*sin(q5)*sin(q6)*cos(q3) - L6*dq5*sin(q2)*sin(q3)*sin(q5)*sin(q6)*cos(q1 + pi/4) + L6*dq5*sin(q2)*sin(q6)*cos(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) + L6*dq5*sin(q4)*sin(q6)*cos(q1 + pi/4)*cos(q2)*cos(q5) + L6*dq6*sin(q1 + pi/4)*sin(q3)*sin(q5)*cos(q4)*cos(q6) - L6*dq6*sin(q1 + pi/4)*cos(q3)*cos(q5)*cos(q6) + L6*dq6*sin(q2)*sin(q3)*cos(q1 + pi/4)*cos(q5)*cos(q6) + L6*dq6*sin(q2)*sin(q5)*cos(q1 + pi/4)*cos(q3)*cos(q4)*cos(q6) + L6*dq6*sin(q4)*sin(q5)*cos(q1 + pi/4)*cos(q2)*cos(q6), 
                       L6*(-dq1*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q5)*cos(q6) - dq1*sin(q1 + pi/4)*sin(q2)*sin(q4)*sin(q6)*cos(q3) + dq1*sin(q1 + pi/4)*sin(q2)*cos(q3)*cos(q4)*cos(q5)*cos(q6) + dq1*sin(q1 + pi/4)*sin(q4)*cos(q2)*cos(q5)*cos(q6) + dq1*sin(q1 + pi/4)*sin(q6)*cos(q2)*cos(q4) + dq1*sin(q3)*sin(q4)*sin(q6)*cos(q1 + pi/4) - dq1*sin(q3)*cos(q1 + pi/4)*cos(q4)*cos(q5)*cos(q6) - dq1*sin(q5)*cos(q1 + pi/4)*cos(q3)*cos(q6) + dq2*sin(q2)*sin(q4)*cos(q1 + pi/4)*cos(q5)*cos(q6) + dq2*sin(q2)*sin(q6)*cos(q1 + pi/4)*cos(q4) + dq2*sin(q3)*sin(q5)*cos(q1 + pi/4)*cos(q2)*cos(q6) + dq2*sin(q4)*sin(q6)*cos(q1 + pi/4)*cos(q2)*cos(q3) - dq2*cos(q1 + pi/4)*cos(q2)*cos(q3)*cos(q4)*cos(q5)*cos(q6) + dq3*sin(q1 + pi/4)*sin(q3)*sin(q5)*cos(q6) + dq3*sin(q1 + pi/4)*sin(q4)*sin(q6)*cos(q3) - dq3*sin(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5)*cos(q6) - dq3*sin(q2)*sin(q3)*sin(q4)*sin(q6)*cos(q1 + pi/4) + dq3*sin(q2)*sin(q3)*cos(q1 + pi/4)*cos(q4)*cos(q5)*cos(q6) + dq3*sin(q2)*sin(q5)*cos(q1 + pi/4)*cos(q3)*cos(q6) + dq4*sin(q1 + pi/4)*sin(q3)*sin(q4)*cos(q5)*cos(q6) + dq4*sin(q1 + pi/4)*sin(q3)*sin(q6)*cos(q4) + dq4*sin(q2)*sin(q4)*cos(q1 + pi/4)*cos(q3)*cos(q5)*cos(q6) + dq4*sin(q2)*sin(q6)*cos(q1 + pi/4)*cos(q3)*cos(q4) + dq4*sin(q4)*sin(q6)*cos(q1 + pi/4)*cos(q2) - dq4*cos(q1 + pi/4)*cos(q2)*cos(q4)*cos(q5)*cos(q6) + dq5*sin(q1 + pi/4)*sin(q3)*sin(q5)*cos(q4)*cos(q6) - dq5*sin(q1 + pi/4)*cos(q3)*cos(q5)*cos(q6) + dq5*sin(q2)*sin(q3)*cos(q1 + pi/4)*cos(q5)*cos(q6) + dq5*sin(q2)*sin(q5)*cos(q1 + pi/4)*cos(q3)*cos(q4)*cos(q6) + dq5*sin(q4)*sin(q5)*cos(q1 + pi/4)*cos(q2)*cos(q6) + dq6*sin(q1 + pi/4)*sin(q3)*sin(q4)*cos(q6) + dq6*sin(q1 + pi/4)*sin(q3)*sin(q6)*cos(q4)*cos(q5) + dq6*sin(q1 + pi/4)*sin(q5)*sin(q6)*cos(q3) - dq6*sin(q2)*sin(q3)*sin(q5)*sin(q6)*cos(q1 + pi/4) + dq6*sin(q2)*sin(q4)*cos(q1 + pi/4)*cos(q3)*cos(q6) + dq6*sin(q2)*sin(q6)*cos(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) + dq6*sin(q4)*sin(q6)*cos(q1 + pi/4)*cos(q2)*cos(q5) - dq6*cos(q1 + pi/4)*cos(q2)*cos(q4)*cos(q6)), 
                       0], 
                      [-L1*dq1*sin(q1 + pi/4) - L2*dq1*sin(q1 + pi/4)*cos(q2) - L2*dq2*sin(q2)*cos(q1 + pi/4) + L3*dq1*sin(q1 + pi/4)*sin(q2)*cos(q3) - L3*dq1*sin(q3)*cos(q1 + pi/4) - L3*dq2*cos(q1 + pi/4)*cos(q2)*cos(q3) - L3*dq3*sin(q1 + pi/4)*cos(q3) + L3*dq3*sin(q2)*sin(q3)*cos(q1 + pi/4) + L4*dq1*sin(q1 + pi/4)*sin(q2)*sin(q4)*cos(q3) - L4*dq1*sin(q1 + pi/4)*cos(q2)*cos(q4) - L4*dq1*sin(q3)*sin(q4)*cos(q1 + pi/4) - L4*dq2*sin(q2)*cos(q1 + pi/4)*cos(q4) - L4*dq2*sin(q4)*cos(q1 + pi/4)*cos(q2)*cos(q3) - L4*dq3*sin(q1 + pi/4)*sin(q4)*cos(q3) + L4*dq3*sin(q2)*sin(q3)*sin(q4)*cos(q1 + pi/4) - L4*dq4*sin(q1 + pi/4)*sin(q3)*cos(q4) - L4*dq4*sin(q2)*cos(q1 + pi/4)*cos(q3)*cos(q4) - L4*dq4*sin(q4)*cos(q1 + pi/4)*cos(q2) - L5*dq1*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q5) + L5*dq1*sin(q1 + pi/4)*sin(q2)*cos(q3)*cos(q4)*cos(q5) + L5*dq1*sin(q1 + pi/4)*sin(q4)*cos(q2)*cos(q5) - L5*dq1*sin(q3)*cos(q1 + pi/4)*cos(q4)*cos(q5) - L5*dq1*sin(q5)*cos(q1 + pi/4)*cos(q3) + L5*dq2*sin(q2)*sin(q4)*cos(q1 + pi/4)*cos(q5) + L5*dq2*sin(q3)*sin(q5)*cos(q1 + pi/4)*cos(q2) - L5*dq2*cos(q1 + pi/4)*cos(q2)*cos(q3)*cos(q4)*cos(q5) + L5*dq3*sin(q1 + pi/4)*sin(q3)*sin(q5) - L5*dq3*sin(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) + L5*dq3*sin(q2)*sin(q3)*cos(q1 + pi/4)*cos(q4)*cos(q5) + L5*dq3*sin(q2)*sin(q5)*cos(q1 + pi/4)*cos(q3) + L5*dq4*sin(q1 + pi/4)*sin(q3)*sin(q4)*cos(q5) + L5*dq4*sin(q2)*sin(q4)*cos(q1 + pi/4)*cos(q3)*cos(q5) - L5*dq4*cos(q1 + pi/4)*cos(q2)*cos(q4)*cos(q5) + L5*dq5*sin(q1 + pi/4)*sin(q3)*sin(q5)*cos(q4) - L5*dq5*sin(q1 + pi/4)*cos(q3)*cos(q5) + L5*dq5*sin(q2)*sin(q3)*cos(q1 + pi/4)*cos(q5) + L5*dq5*sin(q2)*sin(q5)*cos(q1 + pi/4)*cos(q3)*cos(q4) + L5*dq5*sin(q4)*sin(q5)*cos(q1 + pi/4)*cos(q2) - L6*dq1*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q5)*sin(q6) + L6*dq1*sin(q1 + pi/4)*sin(q2)*sin(q4)*cos(q3)*cos(q6) + L6*dq1*sin(q1 + pi/4)*sin(q2)*sin(q6)*cos(q3)*cos(q4)*cos(q5) + L6*dq1*sin(q1 + pi/4)*sin(q4)*sin(q6)*cos(q2)*cos(q5) - L6*dq1*sin(q1 + pi/4)*cos(q2)*cos(q4)*cos(q6) - L6*dq1*sin(q3)*sin(q4)*cos(q1 + pi/4)*cos(q6) - L6*dq1*sin(q3)*sin(q6)*cos(q1 + pi/4)*cos(q4)*cos(q5) - L6*dq1*sin(q5)*sin(q6)*cos(q1 + pi/4)*cos(q3) + L6*dq2*sin(q2)*sin(q4)*sin(q6)*cos(q1 + pi/4)*cos(q5) - L6*dq2*sin(q2)*cos(q1 + pi/4)*cos(q4)*cos(q6) + L6*dq2*sin(q3)*sin(q5)*sin(q6)*cos(q1 + pi/4)*cos(q2) - L6*dq2*sin(q4)*cos(q1 + pi/4)*cos(q2)*cos(q3)*cos(q6) - L6*dq2*sin(q6)*cos(q1 + pi/4)*cos(q2)*cos(q3)*cos(q4)*cos(q5) + L6*dq3*sin(q1 + pi/4)*sin(q3)*sin(q5)*sin(q6) - L6*dq3*sin(q1 + pi/4)*sin(q4)*cos(q3)*cos(q6) - L6*dq3*sin(q1 + pi/4)*sin(q6)*cos(q3)*cos(q4)*cos(q5) + L6*dq3*sin(q2)*sin(q3)*sin(q4)*cos(q1 + pi/4)*cos(q6) + L6*dq3*sin(q2)*sin(q3)*sin(q6)*cos(q1 + pi/4)*cos(q4)*cos(q5) + L6*dq3*sin(q2)*sin(q5)*sin(q6)*cos(q1 + pi/4)*cos(q3) + L6*dq4*sin(q1 + pi/4)*sin(q3)*sin(q4)*sin(q6)*cos(q5) - L6*dq4*sin(q1 + pi/4)*sin(q3)*cos(q4)*cos(q6) + L6*dq4*sin(q2)*sin(q4)*sin(q6)*cos(q1 + pi/4)*cos(q3)*cos(q5) - L6*dq4*sin(q2)*cos(q1 + pi/4)*cos(q3)*cos(q4)*cos(q6) - L6*dq4*sin(q4)*cos(q1 + pi/4)*cos(q2)*cos(q6) - L6*dq4*sin(q6)*cos(q1 + pi/4)*cos(q2)*cos(q4)*cos(q5) + L6*dq5*sin(q1 + pi/4)*sin(q3)*sin(q5)*sin(q6)*cos(q4) - L6*dq5*sin(q1 + pi/4)*sin(q6)*cos(q3)*cos(q5) + L6*dq5*sin(q2)*sin(q3)*sin(q6)*cos(q1 + pi/4)*cos(q5) + L6*dq5*sin(q2)*sin(q5)*sin(q6)*cos(q1 + pi/4)*cos(q3)*cos(q4) + L6*dq5*sin(q4)*sin(q5)*sin(q6)*cos(q1 + pi/4)*cos(q2) + L6*dq6*sin(q1 + pi/4)*sin(q3)*sin(q4)*sin(q6) - L6*dq6*sin(q1 + pi/4)*sin(q3)*cos(q4)*cos(q5)*cos(q6) - L6*dq6*sin(q1 + pi/4)*sin(q5)*cos(q3)*cos(q6) + L6*dq6*sin(q2)*sin(q3)*sin(q5)*cos(q1 + pi/4)*cos(q6) + L6*dq6*sin(q2)*sin(q4)*sin(q6)*cos(q1 + pi/4)*cos(q3) - L6*dq6*sin(q2)*cos(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5)*cos(q6) - L6*dq6*sin(q4)*cos(q1 + pi/4)*cos(q2)*cos(q5)*cos(q6) - L6*dq6*sin(q6)*cos(q1 + pi/4)*cos(q2)*cos(q4), 
                       dq1*(-L2*sin(q2) - L3*cos(q2)*cos(q3) - L4*sin(q2)*cos(q4) - L4*sin(q4)*cos(q2)*cos(q3) + L5*sin(q2)*sin(q4)*cos(q5) + L5*sin(q3)*sin(q5)*cos(q2) - L5*cos(q2)*cos(q3)*cos(q4)*cos(q5) + L6*sin(q2)*sin(q4)*sin(q6)*cos(q5) - L6*sin(q2)*cos(q4)*cos(q6) + L6*sin(q3)*sin(q5)*sin(q6)*cos(q2) - L6*sin(q4)*cos(q2)*cos(q3)*cos(q6) - L6*sin(q6)*cos(q2)*cos(q3)*cos(q4)*cos(q5))*cos(q1 + pi/4) + (-L2*dq2*cos(q2) + L3*dq2*sin(q2)*cos(q3) + L3*dq3*sin(q3)*cos(q2) + L4*dq2*sin(q2)*sin(q4)*cos(q3) - L4*dq2*cos(q2)*cos(q4) + L4*dq3*sin(q3)*sin(q4)*cos(q2) + L4*dq4*sin(q2)*sin(q4) - L4*dq4*cos(q2)*cos(q3)*cos(q4) - L5*dq2*sin(q2)*sin(q3)*sin(q5) + L5*dq2*sin(q2)*cos(q3)*cos(q4)*cos(q5) + L5*dq2*sin(q4)*cos(q2)*cos(q5) + L5*dq3*sin(q3)*cos(q2)*cos(q4)*cos(q5) + L5*dq3*sin(q5)*cos(q2)*cos(q3) + L5*dq4*sin(q2)*cos(q4)*cos(q5) + L5*dq4*sin(q4)*cos(q2)*cos(q3)*cos(q5) - L5*dq5*sin(q2)*sin(q4)*sin(q5) + L5*dq5*sin(q3)*cos(q2)*cos(q5) + L5*dq5*sin(q5)*cos(q2)*cos(q3)*cos(q4) - L6*dq2*sin(q2)*sin(q3)*sin(q5)*sin(q6) + L6*dq2*sin(q2)*sin(q4)*cos(q3)*cos(q6) + L6*dq2*sin(q2)*sin(q6)*cos(q3)*cos(q4)*cos(q5) + L6*dq2*sin(q4)*sin(q6)*cos(q2)*cos(q5) - L6*dq2*cos(q2)*cos(q4)*cos(q6) + L6*dq3*sin(q3)*sin(q4)*cos(q2)*cos(q6) + L6*dq3*sin(q3)*sin(q6)*cos(q2)*cos(q4)*cos(q5) + L6*dq3*sin(q5)*sin(q6)*cos(q2)*cos(q3) + L6*dq4*sin(q2)*sin(q4)*cos(q6) + L6*dq4*sin(q2)*sin(q6)*cos(q4)*cos(q5) + L6*dq4*sin(q4)*sin(q6)*cos(q2)*cos(q3)*cos(q5) - L6*dq4*cos(q2)*cos(q3)*cos(q4)*cos(q6) - L6*dq5*sin(q2)*sin(q4)*sin(q5)*sin(q6) + L6*dq5*sin(q3)*sin(q6)*cos(q2)*cos(q5) + L6*dq5*sin(q5)*sin(q6)*cos(q2)*cos(q3)*cos(q4) + L6*dq6*sin(q2)*sin(q4)*cos(q5)*cos(q6) + L6*dq6*sin(q2)*sin(q6)*cos(q4) + L6*dq6*sin(q3)*sin(q5)*cos(q2)*cos(q6) + L6*dq6*sin(q4)*sin(q6)*cos(q2)*cos(q3) - L6*dq6*cos(q2)*cos(q3)*cos(q4)*cos(q5)*cos(q6))*sin(q1 + pi/4), 
                       -L3*dq1*sin(q1 + pi/4)*cos(q3) + L3*dq1*sin(q2)*sin(q3)*cos(q1 + pi/4) + L3*dq2*sin(q1 + pi/4)*sin(q3)*cos(q2) + L3*dq3*sin(q1 + pi/4)*sin(q2)*cos(q3) - L3*dq3*sin(q3)*cos(q1 + pi/4) - L4*dq1*sin(q1 + pi/4)*sin(q4)*cos(q3) + L4*dq1*sin(q2)*sin(q3)*sin(q4)*cos(q1 + pi/4) + L4*dq2*sin(q1 + pi/4)*sin(q3)*sin(q4)*cos(q2) + L4*dq3*sin(q1 + pi/4)*sin(q2)*sin(q4)*cos(q3) - L4*dq3*sin(q3)*sin(q4)*cos(q1 + pi/4) + L4*dq4*sin(q1 + pi/4)*sin(q2)*sin(q3)*cos(q4) + L4*dq4*cos(q1 + pi/4)*cos(q3)*cos(q4) + L5*dq1*sin(q1 + pi/4)*sin(q3)*sin(q5) - L5*dq1*sin(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5) + L5*dq1*sin(q2)*sin(q3)*cos(q1 + pi/4)*cos(q4)*cos(q5) + L5*dq1*sin(q2)*sin(q5)*cos(q1 + pi/4)*cos(q3) + L5*dq2*sin(q1 + pi/4)*sin(q3)*cos(q2)*cos(q4)*cos(q5) + L5*dq2*sin(q1 + pi/4)*sin(q5)*cos(q2)*cos(q3) - L5*dq3*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q5) + L5*dq3*sin(q1 + pi/4)*sin(q2)*cos(q3)*cos(q4)*cos(q5) - L5*dq3*sin(q3)*cos(q1 + pi/4)*cos(q4)*cos(q5) - L5*dq3*sin(q5)*cos(q1 + pi/4)*cos(q3) - L5*dq4*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q4)*cos(q5) - L5*dq4*sin(q4)*cos(q1 + pi/4)*cos(q3)*cos(q5) - L5*dq5*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q5)*cos(q4) + L5*dq5*sin(q1 + pi/4)*sin(q2)*cos(q3)*cos(q5) - L5*dq5*sin(q3)*cos(q1 + pi/4)*cos(q5) - L5*dq5*sin(q5)*cos(q1 + pi/4)*cos(q3)*cos(q4) + L6*dq1*sin(q1 + pi/4)*sin(q3)*sin(q5)*sin(q6) - L6*dq1*sin(q1 + pi/4)*sin(q4)*cos(q3)*cos(q6) - L6*dq1*sin(q1 + pi/4)*sin(q6)*cos(q3)*cos(q4)*cos(q5) + L6*dq1*sin(q2)*sin(q3)*sin(q4)*cos(q1 + pi/4)*cos(q6) + L6*dq1*sin(q2)*sin(q3)*sin(q6)*cos(q1 + pi/4)*cos(q4)*cos(q5) + L6*dq1*sin(q2)*sin(q5)*sin(q6)*cos(q1 + pi/4)*cos(q3) + L6*dq2*sin(q1 + pi/4)*sin(q3)*sin(q4)*cos(q2)*cos(q6) + L6*dq2*sin(q1 + pi/4)*sin(q3)*sin(q6)*cos(q2)*cos(q4)*cos(q5) + L6*dq2*sin(q1 + pi/4)*sin(q5)*sin(q6)*cos(q2)*cos(q3) - L6*dq3*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q5)*sin(q6) + L6*dq3*sin(q1 + pi/4)*sin(q2)*sin(q4)*cos(q3)*cos(q6) + L6*dq3*sin(q1 + pi/4)*sin(q2)*sin(q6)*cos(q3)*cos(q4)*cos(q5) - L6*dq3*sin(q3)*sin(q4)*cos(q1 + pi/4)*cos(q6) - L6*dq3*sin(q3)*sin(q6)*cos(q1 + pi/4)*cos(q4)*cos(q5) - L6*dq3*sin(q5)*sin(q6)*cos(q1 + pi/4)*cos(q3) - L6*dq4*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q4)*sin(q6)*cos(q5) + L6*dq4*sin(q1 + pi/4)*sin(q2)*sin(q3)*cos(q4)*cos(q6) - L6*dq4*sin(q4)*sin(q6)*cos(q1 + pi/4)*cos(q3)*cos(q5) + L6*dq4*cos(q1 + pi/4)*cos(q3)*cos(q4)*cos(q6) - L6*dq5*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q5)*sin(q6)*cos(q4) + L6*dq5*sin(q1 + pi/4)*sin(q2)*sin(q6)*cos(q3)*cos(q5) - L6*dq5*sin(q3)*sin(q6)*cos(q1 + pi/4)*cos(q5) - L6*dq5*sin(q5)*sin(q6)*cos(q1 + pi/4)*cos(q3)*cos(q4) - L6*dq6*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q4)*sin(q6) + L6*dq6*sin(q1 + pi/4)*sin(q2)*sin(q3)*cos(q4)*cos(q5)*cos(q6) + L6*dq6*sin(q1 + pi/4)*sin(q2)*sin(q5)*cos(q3)*cos(q6) - L6*dq6*sin(q3)*sin(q5)*cos(q1 + pi/4)*cos(q6) - L6*dq6*sin(q4)*sin(q6)*cos(q1 + pi/4)*cos(q3) + L6*dq6*cos(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5)*cos(q6), 
                       -L4*dq1*sin(q1 + pi/4)*sin(q3)*cos(q4) - L4*dq1*sin(q2)*cos(q1 + pi/4)*cos(q3)*cos(q4) - L4*dq1*sin(q4)*cos(q1 + pi/4)*cos(q2) + L4*dq2*sin(q1 + pi/4)*sin(q2)*sin(q4) - L4*dq2*sin(q1 + pi/4)*cos(q2)*cos(q3)*cos(q4) + L4*dq3*sin(q1 + pi/4)*sin(q2)*sin(q3)*cos(q4) + L4*dq3*cos(q1 + pi/4)*cos(q3)*cos(q4) + L4*dq4*sin(q1 + pi/4)*sin(q2)*sin(q4)*cos(q3) - L4*dq4*sin(q1 + pi/4)*cos(q2)*cos(q4) - L4*dq4*sin(q3)*sin(q4)*cos(q1 + pi/4) + L5*dq1*sin(q1 + pi/4)*sin(q3)*sin(q4)*cos(q5) + L5*dq1*sin(q2)*sin(q4)*cos(q1 + pi/4)*cos(q3)*cos(q5) - L5*dq1*cos(q1 + pi/4)*cos(q2)*cos(q4)*cos(q5) + L5*dq2*sin(q1 + pi/4)*sin(q2)*cos(q4)*cos(q5) + L5*dq2*sin(q1 + pi/4)*sin(q4)*cos(q2)*cos(q3)*cos(q5) - L5*dq3*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q4)*cos(q5) - L5*dq3*sin(q4)*cos(q1 + pi/4)*cos(q3)*cos(q5) + L5*dq4*sin(q1 + pi/4)*sin(q2)*cos(q3)*cos(q4)*cos(q5) + L5*dq4*sin(q1 + pi/4)*sin(q4)*cos(q2)*cos(q5) - L5*dq4*sin(q3)*cos(q1 + pi/4)*cos(q4)*cos(q5) - L5*dq5*sin(q1 + pi/4)*sin(q2)*sin(q4)*sin(q5)*cos(q3) + L5*dq5*sin(q1 + pi/4)*sin(q5)*cos(q2)*cos(q4) + L5*dq5*sin(q3)*sin(q4)*sin(q5)*cos(q1 + pi/4) + L6*dq1*sin(q1 + pi/4)*sin(q3)*sin(q4)*sin(q6)*cos(q5) - L6*dq1*sin(q1 + pi/4)*sin(q3)*cos(q4)*cos(q6) + L6*dq1*sin(q2)*sin(q4)*sin(q6)*cos(q1 + pi/4)*cos(q3)*cos(q5) - L6*dq1*sin(q2)*cos(q1 + pi/4)*cos(q3)*cos(q4)*cos(q6) - L6*dq1*sin(q4)*cos(q1 + pi/4)*cos(q2)*cos(q6) - L6*dq1*sin(q6)*cos(q1 + pi/4)*cos(q2)*cos(q4)*cos(q5) + L6*dq2*sin(q1 + pi/4)*sin(q2)*sin(q4)*cos(q6) + L6*dq2*sin(q1 + pi/4)*sin(q2)*sin(q6)*cos(q4)*cos(q5) + L6*dq2*sin(q1 + pi/4)*sin(q4)*sin(q6)*cos(q2)*cos(q3)*cos(q5) - L6*dq2*sin(q1 + pi/4)*cos(q2)*cos(q3)*cos(q4)*cos(q6) - L6*dq3*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q4)*sin(q6)*cos(q5) + L6*dq3*sin(q1 + pi/4)*sin(q2)*sin(q3)*cos(q4)*cos(q6) - L6*dq3*sin(q4)*sin(q6)*cos(q1 + pi/4)*cos(q3)*cos(q5) + L6*dq3*cos(q1 + pi/4)*cos(q3)*cos(q4)*cos(q6) + L6*dq4*sin(q1 + pi/4)*sin(q2)*sin(q4)*cos(q3)*cos(q6) + L6*dq4*sin(q1 + pi/4)*sin(q2)*sin(q6)*cos(q3)*cos(q4)*cos(q5) + L6*dq4*sin(q1 + pi/4)*sin(q4)*sin(q6)*cos(q2)*cos(q5) - L6*dq4*sin(q1 + pi/4)*cos(q2)*cos(q4)*cos(q6) - L6*dq4*sin(q3)*sin(q4)*cos(q1 + pi/4)*cos(q6) - L6*dq4*sin(q3)*sin(q6)*cos(q1 + pi/4)*cos(q4)*cos(q5) - L6*dq5*sin(q1 + pi/4)*sin(q2)*sin(q4)*sin(q5)*sin(q6)*cos(q3) + L6*dq5*sin(q1 + pi/4)*sin(q5)*sin(q6)*cos(q2)*cos(q4) + L6*dq5*sin(q3)*sin(q4)*sin(q5)*sin(q6)*cos(q1 + pi/4) + L6*dq6*sin(q1 + pi/4)*sin(q2)*sin(q4)*cos(q3)*cos(q5)*cos(q6) + L6*dq6*sin(q1 + pi/4)*sin(q2)*sin(q6)*cos(q3)*cos(q4) + L6*dq6*sin(q1 + pi/4)*sin(q4)*sin(q6)*cos(q2) - L6*dq6*sin(q1 + pi/4)*cos(q2)*cos(q4)*cos(q5)*cos(q6) - L6*dq6*sin(q3)*sin(q4)*cos(q1 + pi/4)*cos(q5)*cos(q6) - L6*dq6*sin(q3)*sin(q6)*cos(q1 + pi/4)*cos(q4), 
                       L5*dq1*sin(q1 + pi/4)*sin(q3)*sin(q5)*cos(q4) - L5*dq1*sin(q1 + pi/4)*cos(q3)*cos(q5) + L5*dq1*sin(q2)*sin(q3)*cos(q1 + pi/4)*cos(q5) + L5*dq1*sin(q2)*sin(q5)*cos(q1 + pi/4)*cos(q3)*cos(q4) + L5*dq1*sin(q4)*sin(q5)*cos(q1 + pi/4)*cos(q2) - L5*dq2*sin(q1 + pi/4)*sin(q2)*sin(q4)*sin(q5) + L5*dq2*sin(q1 + pi/4)*sin(q3)*cos(q2)*cos(q5) + L5*dq2*sin(q1 + pi/4)*sin(q5)*cos(q2)*cos(q3)*cos(q4) - L5*dq3*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q5)*cos(q4) + L5*dq3*sin(q1 + pi/4)*sin(q2)*cos(q3)*cos(q5) - L5*dq3*sin(q3)*cos(q1 + pi/4)*cos(q5) - L5*dq3*sin(q5)*cos(q1 + pi/4)*cos(q3)*cos(q4) - L5*dq4*sin(q1 + pi/4)*sin(q2)*sin(q4)*sin(q5)*cos(q3) + L5*dq4*sin(q1 + pi/4)*sin(q5)*cos(q2)*cos(q4) + L5*dq4*sin(q3)*sin(q4)*sin(q5)*cos(q1 + pi/4) - L5*dq5*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q5) + L5*dq5*sin(q1 + pi/4)*sin(q2)*cos(q3)*cos(q4)*cos(q5) + L5*dq5*sin(q1 + pi/4)*sin(q4)*cos(q2)*cos(q5) - L5*dq5*sin(q3)*cos(q1 + pi/4)*cos(q4)*cos(q5) - L5*dq5*sin(q5)*cos(q1 + pi/4)*cos(q3) + L6*dq1*sin(q1 + pi/4)*sin(q3)*sin(q5)*sin(q6)*cos(q4) - L6*dq1*sin(q1 + pi/4)*sin(q6)*cos(q3)*cos(q5) + L6*dq1*sin(q2)*sin(q3)*sin(q6)*cos(q1 + pi/4)*cos(q5) + L6*dq1*sin(q2)*sin(q5)*sin(q6)*cos(q1 + pi/4)*cos(q3)*cos(q4) + L6*dq1*sin(q4)*sin(q5)*sin(q6)*cos(q1 + pi/4)*cos(q2) - L6*dq2*sin(q1 + pi/4)*sin(q2)*sin(q4)*sin(q5)*sin(q6) + L6*dq2*sin(q1 + pi/4)*sin(q3)*sin(q6)*cos(q2)*cos(q5) + L6*dq2*sin(q1 + pi/4)*sin(q5)*sin(q6)*cos(q2)*cos(q3)*cos(q4) - L6*dq3*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q5)*sin(q6)*cos(q4) + L6*dq3*sin(q1 + pi/4)*sin(q2)*sin(q6)*cos(q3)*cos(q5) - L6*dq3*sin(q3)*sin(q6)*cos(q1 + pi/4)*cos(q5) - L6*dq3*sin(q5)*sin(q6)*cos(q1 + pi/4)*cos(q3)*cos(q4) - L6*dq4*sin(q1 + pi/4)*sin(q2)*sin(q4)*sin(q5)*sin(q6)*cos(q3) + L6*dq4*sin(q1 + pi/4)*sin(q5)*sin(q6)*cos(q2)*cos(q4) + L6*dq4*sin(q3)*sin(q4)*sin(q5)*sin(q6)*cos(q1 + pi/4) - L6*dq5*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q5)*sin(q6) + L6*dq5*sin(q1 + pi/4)*sin(q2)*sin(q6)*cos(q3)*cos(q4)*cos(q5) + L6*dq5*sin(q1 + pi/4)*sin(q4)*sin(q6)*cos(q2)*cos(q5) - L6*dq5*sin(q3)*sin(q6)*cos(q1 + pi/4)*cos(q4)*cos(q5) - L6*dq5*sin(q5)*sin(q6)*cos(q1 + pi/4)*cos(q3) + L6*dq6*sin(q1 + pi/4)*sin(q2)*sin(q3)*cos(q5)*cos(q6) + L6*dq6*sin(q1 + pi/4)*sin(q2)*sin(q5)*cos(q3)*cos(q4)*cos(q6) + L6*dq6*sin(q1 + pi/4)*sin(q4)*sin(q5)*cos(q2)*cos(q6) - L6*dq6*sin(q3)*sin(q5)*cos(q1 + pi/4)*cos(q4)*cos(q6) + L6*dq6*cos(q1 + pi/4)*cos(q3)*cos(q5)*cos(q6), 
                       L6*(dq1*sin(q1 + pi/4)*sin(q3)*sin(q4)*sin(q6) - dq1*sin(q1 + pi/4)*sin(q3)*cos(q4)*cos(q5)*cos(q6) - dq1*sin(q1 + pi/4)*sin(q5)*cos(q3)*cos(q6) + dq1*sin(q2)*sin(q3)*sin(q5)*cos(q1 + pi/4)*cos(q6) + dq1*sin(q2)*sin(q4)*sin(q6)*cos(q1 + pi/4)*cos(q3) - dq1*sin(q2)*cos(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5)*cos(q6) - dq1*sin(q4)*cos(q1 + pi/4)*cos(q2)*cos(q5)*cos(q6) - dq1*sin(q6)*cos(q1 + pi/4)*cos(q2)*cos(q4) + dq2*sin(q1 + pi/4)*sin(q2)*sin(q4)*cos(q5)*cos(q6) + dq2*sin(q1 + pi/4)*sin(q2)*sin(q6)*cos(q4) + dq2*sin(q1 + pi/4)*sin(q3)*sin(q5)*cos(q2)*cos(q6) + dq2*sin(q1 + pi/4)*sin(q4)*sin(q6)*cos(q2)*cos(q3) - dq2*sin(q1 + pi/4)*cos(q2)*cos(q3)*cos(q4)*cos(q5)*cos(q6) - dq3*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q4)*sin(q6) + dq3*sin(q1 + pi/4)*sin(q2)*sin(q3)*cos(q4)*cos(q5)*cos(q6) + dq3*sin(q1 + pi/4)*sin(q2)*sin(q5)*cos(q3)*cos(q6) - dq3*sin(q3)*sin(q5)*cos(q1 + pi/4)*cos(q6) - dq3*sin(q4)*sin(q6)*cos(q1 + pi/4)*cos(q3) + dq3*cos(q1 + pi/4)*cos(q3)*cos(q4)*cos(q5)*cos(q6) + dq4*sin(q1 + pi/4)*sin(q2)*sin(q4)*cos(q3)*cos(q5)*cos(q6) + dq4*sin(q1 + pi/4)*sin(q2)*sin(q6)*cos(q3)*cos(q4) + dq4*sin(q1 + pi/4)*sin(q4)*sin(q6)*cos(q2) - dq4*sin(q1 + pi/4)*cos(q2)*cos(q4)*cos(q5)*cos(q6) - dq4*sin(q3)*sin(q4)*cos(q1 + pi/4)*cos(q5)*cos(q6) - dq4*sin(q3)*sin(q6)*cos(q1 + pi/4)*cos(q4) + dq5*sin(q1 + pi/4)*sin(q2)*sin(q3)*cos(q5)*cos(q6) + dq5*sin(q1 + pi/4)*sin(q2)*sin(q5)*cos(q3)*cos(q4)*cos(q6) + dq5*sin(q1 + pi/4)*sin(q4)*sin(q5)*cos(q2)*cos(q6) - dq5*sin(q3)*sin(q5)*cos(q1 + pi/4)*cos(q4)*cos(q6) + dq5*cos(q1 + pi/4)*cos(q3)*cos(q5)*cos(q6) - dq6*sin(q1 + pi/4)*sin(q2)*sin(q3)*sin(q5)*sin(q6) + dq6*sin(q1 + pi/4)*sin(q2)*sin(q4)*cos(q3)*cos(q6) + dq6*sin(q1 + pi/4)*sin(q2)*sin(q6)*cos(q3)*cos(q4)*cos(q5) + dq6*sin(q1 + pi/4)*sin(q4)*sin(q6)*cos(q2)*cos(q5) - dq6*sin(q1 + pi/4)*cos(q2)*cos(q4)*cos(q6) - dq6*sin(q3)*sin(q4)*cos(q1 + pi/4)*cos(q6) - dq6*sin(q3)*sin(q6)*cos(q1 + pi/4)*cos(q4)*cos(q5) - dq6*sin(q5)*sin(q6)*cos(q1 + pi/4)*cos(q3)), 0], 
                      [0, 
                       L2*dq2*sin(q2) + L3*dq2*cos(q2)*cos(q3) - L3*dq3*sin(q2)*sin(q3) + L4*(dq2*sin(q2)*cos(q4) + dq2*sin(q4)*cos(q2)*cos(q3) - dq3*sin(q2)*sin(q3)*sin(q4) + dq4*sin(q2)*cos(q3)*cos(q4) + dq4*sin(q4)*cos(q2)) + L5*(-dq2*sin(q3)*sin(q5)*cos(q2) - dq3*sin(q2)*sin(q5)*cos(q3) - dq5*(sin(q2)*cos(q3)*cos(q4) + sin(q4)*cos(q2))*sin(q5) - dq5*sin(q2)*sin(q3)*cos(q5) + (-dq2*sin(q2)*sin(q4) + dq2*cos(q2)*cos(q3)*cos(q4) - dq3*sin(q2)*sin(q3)*cos(q4) - dq4*sin(q2)*sin(q4)*cos(q3) + dq4*cos(q2)*cos(q4))*cos(q5)) + L6*(dq6*((sin(q2)*cos(q3)*cos(q4) + sin(q4)*cos(q2))*cos(q5) - sin(q2)*sin(q3)*sin(q5))*cos(q6) - dq6*(sin(q2)*sin(q4)*cos(q3) - cos(q2)*cos(q4))*sin(q6) + (dq2*sin(q2)*cos(q4) + dq2*sin(q4)*cos(q2)*cos(q3) - dq3*sin(q2)*sin(q3)*sin(q4) + dq4*sin(q2)*cos(q3)*cos(q4) + dq4*sin(q4)*cos(q2))*cos(q6) + (-dq2*sin(q3)*sin(q5)*cos(q2) - dq3*sin(q2)*sin(q5)*cos(q3) - dq5*(sin(q2)*cos(q3)*cos(q4) + sin(q4)*cos(q2))*sin(q5) - dq5*sin(q2)*sin(q3)*cos(q5) + (-dq2*sin(q2)*sin(q4) + dq2*cos(q2)*cos(q3)*cos(q4) - dq3*sin(q2)*sin(q3)*cos(q4) - dq4*sin(q2)*sin(q4)*cos(q3) + dq4*cos(q2)*cos(q4))*cos(q5))*sin(q6)), 
                       -dq2*(L3*sin(q3) + L4*sin(q3)*sin(q4) + L5*(sin(q3)*cos(q4)*cos(q5) + sin(q5)*cos(q3)) + L6*((sin(q3)*cos(q4)*cos(q5) + sin(q5)*cos(q3))*sin(q6) + sin(q3)*sin(q4)*cos(q6)))*sin(q2) + (L3*dq3*cos(q3) + L4*dq3*sin(q4)*cos(q3) + L4*dq4*sin(q3)*cos(q4) + L5*(-dq3*sin(q3)*sin(q5) + dq3*cos(q3)*cos(q4)*cos(q5) - dq4*sin(q3)*sin(q4)*cos(q5) - dq5*sin(q3)*sin(q5)*cos(q4) + dq5*cos(q3)*cos(q5)) + L6*(dq3*sin(q4)*cos(q3)*cos(q6) + dq4*sin(q3)*cos(q4)*cos(q6) + dq6*(sin(q3)*cos(q4)*cos(q5) + sin(q5)*cos(q3))*cos(q6) - dq6*sin(q3)*sin(q4)*sin(q6) + (-dq3*sin(q3)*sin(q5) + dq3*cos(q3)*cos(q4)*cos(q5) - dq4*sin(q3)*sin(q4)*cos(q5) - dq5*sin(q3)*sin(q5)*cos(q4) + dq5*cos(q3)*cos(q5))*sin(q6)))*cos(q2), 
                       L4*(dq2*sin(q2)*cos(q3)*cos(q4) + dq2*sin(q4)*cos(q2) + dq3*sin(q3)*cos(q2)*cos(q4) + dq4*sin(q2)*cos(q4) + dq4*sin(q4)*cos(q2)*cos(q3)) - L5*dq5*(sin(q2)*cos(q4) + sin(q4)*cos(q2)*cos(q3))*sin(q5) + L5*(-dq2*sin(q2)*sin(q4)*cos(q3) + dq2*cos(q2)*cos(q4) - dq3*sin(q3)*sin(q4)*cos(q2) - dq4*sin(q2)*sin(q4) + dq4*cos(q2)*cos(q3)*cos(q4))*cos(q5) + L6*(-dq5*(sin(q2)*cos(q4) + sin(q4)*cos(q2)*cos(q3))*sin(q5)*sin(q6) - dq6*(sin(q2)*sin(q4) - cos(q2)*cos(q3)*cos(q4))*sin(q6) + dq6*(sin(q2)*cos(q4) + sin(q4)*cos(q2)*cos(q3))*cos(q5)*cos(q6) + (-dq2*sin(q2)*sin(q4)*cos(q3) + dq2*cos(q2)*cos(q4) - dq3*sin(q3)*sin(q4)*cos(q2) - dq4*sin(q2)*sin(q4) + dq4*cos(q2)*cos(q3)*cos(q4))*sin(q6)*cos(q5) + (dq2*sin(q2)*cos(q3)*cos(q4) + dq2*sin(q4)*cos(q2) + dq3*sin(q3)*cos(q2)*cos(q4) + dq4*sin(q2)*cos(q4) + dq4*sin(q4)*cos(q2)*cos(q3))*cos(q6)), 
                       -L6*dq6*((sin(q2)*sin(q4) - cos(q2)*cos(q3)*cos(q4))*sin(q5) - sin(q3)*cos(q2)*cos(q5))*cos(q6) + (-L5 - L6*sin(q6))*(dq2*sin(q2)*sin(q3)*cos(q5) - dq3*cos(q2)*cos(q3)*cos(q5) + dq5*(sin(q2)*sin(q4) - cos(q2)*cos(q3)*cos(q4))*cos(q5) + dq5*sin(q3)*sin(q5)*cos(q2) + (dq2*sin(q2)*cos(q3)*cos(q4) + dq2*sin(q4)*cos(q2) + dq3*sin(q3)*cos(q2)*cos(q4) + dq4*sin(q2)*cos(q4) + dq4*sin(q4)*cos(q2)*cos(q3))*sin(q5)), 
                       L6*(-dq6*((sin(q2)*sin(q4) - cos(q2)*cos(q3)*cos(q4))*cos(q5) + sin(q3)*sin(q5)*cos(q2))*sin(q6) + dq6*(sin(q2)*cos(q4) + sin(q4)*cos(q2)*cos(q3))*cos(q6) + (-dq2*sin(q2)*sin(q3)*sin(q5) + dq3*sin(q5)*cos(q2)*cos(q3) - dq5*(sin(q2)*sin(q4) - cos(q2)*cos(q3)*cos(q4))*sin(q5) + dq5*sin(q3)*cos(q2)*cos(q5) + (dq2*sin(q2)*cos(q3)*cos(q4) + dq2*sin(q4)*cos(q2) + dq3*sin(q3)*cos(q2)*cos(q4) + dq4*sin(q2)*cos(q4) + dq4*sin(q4)*cos(q2)*cos(q3))*cos(q5))*cos(q6) + (-dq2*sin(q2)*sin(q4)*cos(q3) + dq2*cos(q2)*cos(q4) - dq3*sin(q3)*sin(q4)*cos(q2) - dq4*sin(q2)*sin(q4) + dq4*cos(q2)*cos(q3)*cos(q4))*sin(q6)), 
                       0]]) 
        return z
    
    def djacobi_all(self, q, dq):
        """ヤコビ行列の時間微分を全部計算する
        ・
        """
        djWo = self.djacobi_Wo(q, dq)
        djBL = self.djacobi_BL(q, dq)
        dj0 = self.djacobi_0(q, dq)
        dj1 = self.djacobi_1(q, dq)
        dj2 = self.djacobi_2(q, dq)
        dj3 = self.djacobi_3(q, dq)
        dj4 = self.djacobi_4(q, dq)
        dj5 = self.djacobi_5(q, dq)
        dj6 = self.djacobi_6(q, dq)
        dj7 = self.djacobi_7(q, dq)
        djGL = self.djacobi_GL(q, dq)
        return [djWo, djBL, dj0, dj1, dj2, dj3, dj4, dj5, dj6, dj7, djGL]

