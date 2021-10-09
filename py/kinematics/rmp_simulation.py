"""RMPシミュレーション"""

#from collections import deque
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anm
import scipy.integrate as integrate

import tqdm
import time

import environment
from new import BaxterRobotArmKinematics
from rmp import RMPfromGDS

from rmp import OriginalRMP


class PositinData:
    def __init__(self,):
        self.x = []
        self.y = []
        self.z = []
    
    def add(self, X):
        self.x.append(X[0, 0])
        self.y.append(X[1, 0])
        self.z.append(X[2, 0])


class FrameData:
    def __init__(self, arm):
        self._get_joint_position(arm)
        self._get_cpoints_position(arm)
    
    
    def _get_joint_position(self, arm):
        self.joint_positions_list = PositinData()
        for o in arm.get_joint_positions():
            self.joint_positions_list.add(o)
        return
    
    
    def _get_cpoints_position(self, arm):
        self.cpoints_potisions_list = []
        for cpoints in arm.cpoints_x:
            _cs = PositinData()
            for c in cpoints:
                _cs.add(c)
            self.cpoints_potisions_list.append(_cs)
        return



class SimulationData:
    
    def __init__(self,):
        self.data = []
        self.ee = PositinData()
        self.command = []
        pass
    
    
    def add_data(self, arm, ddq=None):
        self.data.append(FrameData(arm))
        
        self.ee.add(arm.Ts_Wo[-1].o)
        if ddq is not None:
            self.command.append(ddq)
        return


def _calc_goal_position(t):
    omega = 0.3
    r = 0.5
    g = np.array([[
        r * np.cos(omega*t),
        0,
        r * np.sin(omega * t),
    ]]).T
    g0 = np.array([[0.3, -0.6, 1]]).T
    return g + g0


def _calc_obstacle_position(t):
    return environment.set_obstacle(data=environment.data1)




class Simulator:
    """"""
    
    def __init__(self, isleft=True, TIME_SPAN=30, TIME_INTERVAL=0.05):
        self.isleft = isleft
        self.TIME_SPAN = TIME_SPAN
        self.TIME_INTERVAL = TIME_INTERVAL
        
        return
    
    
    def set_controller(self, rmp_param):
        """rmpをセット"""
        
        self.rmps = []
        for i in range(7):
            if rmp_param[i]['name'] == 'origimal':
                self.rmps.append(OriginalRMP(**rmp_param[i]))
            elif rmp_param[i]['name'] == 'fromGDS':
                self.rmps.append(RMPfromGDS(**rmp_param[i]))
    
        return
    
    
    
    def set_environment(self, env_param):
        
        
        
        
        
        return
    
    
    def run_simulation(self,):
        
        
        self.gl_goal = np.array([[0.3, -0.75, 1]]).T
        #self.obs = [np.array([[0.8, -0.6, 1]]).T]
        
        
        # self.obs = [
        #     np.array([[0.6, -0.6, 1]]).T,
        #     np.array([[0.6, -0.6, 1.1]]).T,
        #     np.array([[0.6, -0.6, 0.9]]).T,
        #     np.array([[0.6, -0.6, 1.2]]).T,
        # ]
        
        
        # self.obs = environment.make_obstacle(name='curb', R=0.15, center=np.array([[0.6, -0.6, 1]]).T)
        # self.obs.extend(make_obstacle(name='curb', R=0.15, center=np.array([[0.6, -0.6, 1.5]]).T))
        # self.obs_plot = np.concatenate(self.obs, axis=1)
        
        
        self.obs = environment.set_obstacle(data=environment.data2)
        
        dobs = np.zeros((3, 1))
        
        #self.obs = None
        
        if self.obs is not None:
            self.obs_plot = np.concatenate(self.obs, axis=1)
        
        t = np.arange(0.0, self.TIME_SPAN, self.TIME_INTERVAL)
        print("t size = ", t.shape)
        
        arm = BaxterRobotArmKinematics(self.isLeft)
        rmp = OriginalRMP(
            attract_max_speed = 2, 
            attract_gain = 10, 
            attract_a_damp_r = 0.15,
            attract_sigma_W = 1, 
            attract_sigma_H = 1, 
            attract_A_damp_r = 5, 
            obs_scale_rep = 0.2,
            obs_scale_damp = 1,
            obs_ratio = 0.5,
            obs_rep_gain = 0.5e-10,
            obs_r = 1,
            jl_gamma_p = 0.05,
            jl_gamma_d = 0.1,
            jl_lambda = 0.7,
            joint_limit_upper = arm.q_max,
            joint_limit_lower = arm.q_min,
        )
        
        rmp2 = RMPfromGDS(
            attract_max_speed = 2, 
            attract_gain = 100,
            attract_alpha_f = 0.3,
            attract_sigma_alpha = 1,
            attract_sigma_gamma = 1,
            attract_w_u = 1,
            attract_w_l = 1,
            attract_alpha = 0.1,
            attract_epsilon = 1e-5,
            jl_gamma_p = 0.05,
            jl_gamma_d = 0.1,
            jl_lambda = 0.7,
            joint_limit_upper = arm.q_max,
            joint_limit_lower = arm.q_min,
            jl_sigma = 1,
        )
        
        
        
        #self.data = SimulationData()
        
        
        def _eom(t, state):
            
            if t > 1 and (int(t) % 2 == 0):
                print("t = ", '{:.2f}'.format(t))
            
            
            #self.gl_goal = _moveing_goal(t)
            q = np.array([state[0:7]]).T
            dq = np.array([state[7:14]]).T
            
            # print("q = ", q)
            # print("dq = ", dq)
            
            arm.update_all(q, dq)  # ロボットアームの全情報更新
            
            pulled_f_all = []
            pulled_M_all = []
            
            for i in range(8):
                for x, dx, J, dJ, rmp in zip(
                    arm.cpoints_x[i],
                    arm.cpoints_dx[i],
                    arm.Jos_cpoints[i],
                    arm.Jos_cpoints_diff_by_t[i],
                    self.rmps[i],
                ):
                    
                    if self.obs is not None:
                        for o in self.obs:
                            a = rmp.a_obs(x, dx, o)
                            M = rmp.metric_obs(x, dx, o, a)
                            f = M @ a

                            
                            
                            _pulled_f = J.T @ (f - M @ dJ @ dq)
                            _pulled_M = J.T @ M @ J
                            
                            pulled_f_all.append(_pulled_f)
                            pulled_M_all.append(_pulled_M)

                    if i == 7:
                        a = rmp.a_attract(x, dx, self.gl_goal)
                        M = rmp.metric_attract(x, dx, self.gl_goal, a)
                        f = M @ a
                        
                        
                        #M = rmp2.inertia_attract(x, dx, self.gl_goal, np.zeros((3,1)))
                        #f = rmp2.f_attract(x, dx, self.gl_goal, np.zeros((3,1)),M)
                        
                        
                        
                        _pulled_f = J.T @ (f - M @ dJ @ dq)
                        _pulled_M = J.T @ M @ J
                        
                        pulled_f_all.append(_pulled_f)
                        pulled_M_all.append(_pulled_M)
            
            pulled_f_all = np.sum(pulled_f_all, axis=0)
            pulled_M_all = np.sum(pulled_M_all, axis=0)
            
            # # ジョイント制限
            # a_jl = rmp.a_joint_limit(q, dq)
            # M_jl = rmp.metric_joint_limit(q)
            # f_jl = M_jl @ a_jl
            # pulled_f_all += f_jl
            # pulled_M_all += M_jl
            
            
            ddq = np.linalg.pinv(pulled_M_all) @ pulled_f_all
            
            dstate = np.concatenate([dq, ddq], axis=0)
            dstate = np.ravel(dstate).tolist()
            
            
            # 以下 視覚化のためのデータ保存
            #self.data.add_data(arm, ddq)
            
            return dstate
        
        
        print("シミュレーション実行中...")
        start = time.time()
        
        # scipy使用
        self.sol = integrate.solve_ivp(
            fun=_eom,
            t_span=(0.0, self.TIME_SPAN),
            y0=np.ravel(np.concatenate([arm.q, arm.dq])).tolist(),
            method='RK45',
            #method='LSODA',
            t_eval=t,
        )
        print("シミュレーション実行終了")
        print("シミュレーション実行時間 = ", time.time() - start)
        print("")
        
        # データ作成
        print("データ作成中...")
        start = time.time()
        self.data = SimulationData()
        _arm = BaxterRobotArmKinematics(isLeft=True)
        
        for i in tqdm.tqdm(range(len(self.sol.t))):
        #for i in range(len(self.sol.t)):
            q = np.array([
                [self.sol.y[0][i]],
                [self.sol.y[1][i]],
                [self.sol.y[2][i]],
                [self.sol.y[3][i]],
                [self.sol.y[4][i]],
                [self.sol.y[5][i]],
                [self.sol.y[6][i]],
            ])
            dq = np.array([
                [self.sol.y[7][i]],
                [self.sol.y[8][i]],
                [self.sol.y[9][i]],
                [self.sol.y[10][i]],
                [self.sol.y[11][i]],
                [self.sol.y[12][i]],
                [self.sol.y[13][i]],
            ])
            _arm.update_all(q, dq)
            self.data.add_data(_arm, ddq=None)
        
        print("データ作成完了")
        print("データ作成時間 = ", time.time() - start)
        print("")
        
        # # オイラー
        # state = np.ravel(np.concatenate([arm.q, arm.dq])).tolist()
        # for i in t:
        #     dstate = _eom(i, state)
        #     for j in range(14):
        #         state[j] = state[j] + dstate[j] * self.TIME_INTERVAL
        


        return




    def plot_animation_2(self,):
        """グラフ作成（遅いかも）"""
        
        start = time.time()
        print("plot実行中...")
        
        
        
        
        # アニメーション
        fig_ani = plt.figure()
        ax = fig_ani.add_subplot(projection = '3d')
        ax.grid(True)
        ax.set_xlabel('X[m]')
        ax.set_ylabel('Y[m]')
        ax.set_zlabel('Z[m]')

        ## 三軸のスケールを揃える
        max_x = 1.0
        min_x = -1.0
        max_y = 0.2
        min_y = -1.0
        max_z = 2.0
        min_z = 0.0
        
        max_range = np.array([
            max_x - min_x,
            max_y - min_y,
            max_z - min_z
            ]).max() * 0.5
        mid_x = (max_x + min_x) * 0.5
        mid_y = (max_y + min_y) * 0.5
        mid_z = (max_z + min_z) * 0.5
        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)


        # 時刻表示
        time_template = 'time = %s [s]'
        

        ax.set_box_aspect((1,1,1))

        def _update(i):
            """アニメーションの関数"""
            ax.cla()  # 遅いかも
            ax.grid(True)
            ax.set_xlabel('X[m]')
            ax.set_ylabel('Y[m]')
            ax.set_zlabel('Z[m]')
            
            ax.set_xlim(mid_x - max_range, mid_x + max_range)
            ax.set_ylim(mid_y - max_range, mid_y + max_range)
            ax.set_zlim(mid_z - max_range, mid_z + max_range)
            
            # 目標点
            ax.scatter(
                self.gl_goal[0, 0], self.gl_goal[1, 0], self.gl_goal[2, 0],
                s = 100, label = 'goal point', marker = '*', color = '#ff7f00', 
                alpha = 1, linewidths = 1.5, edgecolors = 'red')
            
            # _g = _moveing_goal(i*self.TIME_INTERVAL)
            # ax.scatter(
            #     _g[0, 0], _g[1, 0], _g[2, 0],
            #     s = 100, label = 'goal point', marker = '*', color = '#ff7f00', 
            #     alpha = 1, linewidths = 1.5, edgecolors = 'red')
            
            
            
            # 障害物点
            if self.obs is not None:
                ax.scatter(
                    self.obs_plot[0, :], self.obs_plot[1, :], self.obs_plot[2, :],
                    label = 'obstacle point', marker = '.', color = 'k',)
            
            d = self.data.data[i]
            
            # ジョイント位置
            ax.plot(
                d.joint_positions_list.x,
                d.joint_positions_list.y,
                d.joint_positions_list.z,
                "o-", color = "blue",
            )
            
            # 制御点
            for p in d.cpoints_potisions_list:
                ax.scatter(
                    p.x, p.y, p.z,
                    marker='o'
                )
            
            # グリッパー
            ax.plot(
                self.data.ee.x[0:i],
                self.data.ee.y[0:i],
                self.data.ee.z[0:i],
                "-", color = "#ff7f00"
            )
            
            # 時刻表示
            ax.text(
                0.8, 0.12, 0.01,
                time_template % (i * self.TIME_INTERVAL), size = 10
            )
            
            ax.set_box_aspect((1,1,1))

            return

        ani = anm.FuncAnimation(
            fig = fig_ani,
            func = _update,
            frames = len(self.data.data),
            #frames = int(self.TIME_SPAN / self.TIME_INTERVAL),
            interval = self.TIME_INTERVAL * 0.001
        )

        #ani.save("hoge.gif", fps=1/self.TIME_INTERVAL, writer='pillow')
        
        
        print("plot実行終了")
        print("実行時間 = ", time.time() - start)
        
        
        
        # # 入力履歴
        # if len(self.data.command) != 1:
        #     fig_input = plt.figure()
        #     ax2 = fig_input.add_subplot(111)
        #     _c = np.concatenate(self.data.command, axis=1)
        #     for i in range(7):
        #         ax2.plot(_c[i, :], label=str(i+1))
        #     ax2.grid(True)
        #     ax2.legend()
        
        
        # 最終結果
        fig_rezult = plt.figure()
        ax_rezult = fig_rezult.add_subplot(projection='3d')
        ax_rezult.grid(True)
        ax_rezult.set_xlabel('X[m]')
        ax_rezult.set_ylabel('Y[m]')
        ax_rezult.set_zlabel('Z[m]')

        ## 三軸のスケールを揃える
        max_x = 1.0
        min_x = -1.0
        max_y = 0.2
        min_y = -1.0
        max_z = 2.0
        min_z = 0.0
        
        max_range = np.array([
            max_x - min_x,
            max_y - min_y,
            max_z - min_z
            ]).max() * 0.5
        mid_x = (max_x + min_x) * 0.5
        mid_y = (max_y + min_y) * 0.5
        mid_z = (max_z + min_z) * 0.5
        ax_rezult.set_xlim(mid_x - max_range, mid_x + max_range)
        ax_rezult.set_ylim(mid_y - max_range, mid_y + max_range)
        ax_rezult.set_zlim(mid_z - max_range, mid_z + max_range)
        
        i = int(self.TIME_SPAN/self.TIME_INTERVAL)-1
        
        #目標点
        ax_rezult.scatter(
            self.gl_goal[0, 0], self.gl_goal[1, 0], self.gl_goal[2, 0],
            s = 100, label = 'goal point', marker = '*', color = '#ff7f00', 
            alpha = 1, linewidths = 1.5, edgecolors = 'red')
        
        # _g = _moveing_goal(i*self.TIME_INTERVAL)
        # ax.scatter(
        #     _g[0, 0], _g[1, 0], _g[2, 0],
        #     s = 100, label = 'goal point', marker = '*', color = '#ff7f00', 
        #     alpha = 1, linewidths = 1.5, edgecolors = 'red')
        
        
        
        # 障害物点
        if self.obs is not None:
            ax_rezult.scatter(
                self.obs_plot[0, :], self.obs_plot[1, :], self.obs_plot[2, :],
                label = 'obstacle point', marker = '.', color = 'k',)
        
        d = self.data.data[i]
        
        # ジョイント位置
        ax_rezult.plot(
            d.joint_positions_list.x,
            d.joint_positions_list.y,
            d.joint_positions_list.z,
            "o-", color = "blue",
        )
        
        # 制御点
        for p in d.cpoints_potisions_list:
            ax_rezult.scatter(
                p.x, p.y, p.z,
                marker='o'
            )
        
        # グリッパー
        ax_rezult.plot(
            self.data.ee.x[0:i],
            self.data.ee.y[0:i],
            self.data.ee.z[0:i],
            "-", color = "#ff7f00"
        )
        
        # 時刻表示
        ax_rezult.text(
            0.8, 0.12, 0.01,
            time_template % (i * self.TIME_INTERVAL), size = 10
        )
        
        
        
        ax_rezult.set_box_aspect((1,1,1))
        
        
        
        plt.show()
        
        return




def main():
    simu = Simulator()
    simu.run_simulation()
    simu.plot_animation_2()
    
    




if __name__ == "__main__":
    main()
