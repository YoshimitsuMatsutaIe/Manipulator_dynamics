"""RMPシミュレーション"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate


import time


from new import BaxterRobotArmKinematics
from rmp import OriginalRMP


class Simulator:
    """"""
    
    def __init__(self,):
        return
    



class Data:
    
    def __init__(self, t, x):
        self.t = t
        self.x = x



def simulation():
    
    
    gl_goal = np.array([[0, 0, 0]]).T
    obs = np.array([[10, 10, 10]]).T
    dobs = np.zeros((3, 1))
    
    TIME_INTERVAL = 0.01
    TIME_SPAN = 10
    
    t = np.arange(0.0, TIME_SPAN, TIME_INTERVAL)
    
    arm = BaxterRobotArmKinematics(isLeft=True)
    rmp = OriginalRMP(
        attract_max_speed = 2, 
        attract_gain = 10, 
        attract_a_damp_r = 0.3,
        attract_sigma_W = 1, 
        attract_sigma_H = 1, 
        attract_A_damp_r = 0.3, 
        obs_scale_rep = 0.2,
        obs_scale_damp = 1,
        obs_ratio = 0.5,
        obs_rep_gain = 0.5,
        obs_r = 15,
        jl_gamma_p = 0.05,
        jl_gamma_d = 0.1,
        jl_lambda = 0.7,
        joint_limit_upper = arm.q_max,
        joint_limit_lower = arm.q_min,
    )
    
    def eom(t, state):
        
        q = np.array([state[0:7]]).T
        dq = np.array([state[7:14]]).T
        
        arm.update_all(q, dq)  # ロボットアームの全情報更新
        
        pulled_f_all = []
        pulled_M_all = []
        
        for i in range(8):
            for x, dx, J in zip(arm.cpoints_x[i], arm.cpoints_dx[i], arm.Jos_cpoints[i]):
                a = rmp.a_obs(x, dx, obs)
                M = rmp.metric_obs(x, dx, obs, a)
                f = M @ a
                
                _pulled_f = J.T @ f
                _pulled_M = J.T @ M @ J
                
                pulled_f_all.append(_pulled_f)
                pulled_M_all.append(_pulled_M)

                if i == 7:
                    a = rmp.a_attract(x, dx, gl_goal)
                    M = rmp.metric_attract(x, dx, gl_goal, a)
                    f = M @ a
                    
                    _pulled_f = J.T @ f
                    _pulled_M = J.T @ M @ J
                    
                    pulled_f_all.append(_pulled_f)
                    pulled_M_all.append(_pulled_M)
        
        pulled_f_all = np.sum(pulled_f_all, axis=0)
        pulled_M_all = np.sum(pulled_M_all, axis=0)
        
        ddq = np.linalg.pinv(pulled_M_all) @ pulled_f_all
        
        dstate = np.concatenate([dq, ddq], axis=0)
        return np.ravel(dstate).tolist()
    
    print("シミュレーション実行中...")
    start = time.time()
    sol = integrate.solve_ivp(
        fun=eom,
        t_span=(0.0, TIME_SPAN),
        y0=np.ravel(np.concatenate([arm.q, arm.dq])).tolist(),
        method='RK45',
        t_eval=t,
    )
    print("シミュレーション実行終了")
    print("実行時間 = ", time.time() - start)
    
    return sol


def main():
    sol = simulation()
    
    
    fig = plt.figure()
    ax = fig.add_subplot()
    for i in range(7):
        ax.plot(sol.t, sol.y[i], label=str(i+1))
    ax.legend()
    ax.grid(True)
    ax.set_xlabel('time')
    ax.set_ylabel('joint angle')
    
    
    arm = BaxterRobotArmKinematics(isLeft=True)
    for i in range(len(sol.t)):
        _q = [sol.y[j][i] for j in range(7)]
        q = np.array([_q]).T
        arm.update_all(q, arm.dq)
        
        
    
    
    
    plt.show()

if __name__ == "__main__":
    print("実行中...")
    start = time.time()
    
    main()
    
    print("実行終了")
    print("実行時間 = ", time.time() - start)