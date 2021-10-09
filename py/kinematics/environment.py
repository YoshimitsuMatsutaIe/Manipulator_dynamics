"""シミュレーション環境"""


import numpy as np
from math import cos, sin, tan, pi
import matplotlib.pyplot as plt

from new import BaxterRobotArmKinematics


def _rotate(theta, phi, zeta):
    """3次元回転行列"""
    
    Rx = np.array([
        [1, 0, 0],
        [0, cos(theta), -sin(theta)],
        [0, sin(theta), cos(theta)],
    ])
    Ry = np.array([
        [cos(phi), 0, sin(phi)],
        [0, 1, 0],
        [-sin(phi), 0, cos(phi)],
    ])
    Rz = np.array([
        [cos(zeta), -sin(zeta), 0],
        [sin(zeta), cos(zeta), 0],
        [0, 0, 1],
    ])
    
    return Rx @ Ry @ Rz


def _set_sphere(r, center, n):
    """
    
    r : 半径
    center : 中心
    n : 点の数
    """
    
    obs = []
    rand = np.random.RandomState(123)
    for i in range(n):
        theta = np.arccos(rand.uniform(-1, 1))
        phi = rand.uniform(0, 2*pi)
        x = r * sin(theta) * cos(phi)
        y = r * sin(theta) * sin(phi)
        z = r * cos(theta)
        obs.append(np.array([[x, y, z]]).T + center)
    
    return obs


def _set_cylinder(r, L, center, n, theta=0, phi=0, zeta=0,):
    """円筒を設置
    
    r : 半径
    L : 長さ
    n : 点の数
    theta : 回転
    phi : 回転
    zeta : 回転
    """
    
    obs = []
    
    rand = np.random.RandomState(123)
    for i in range(n):
        _theta = rand.uniform(0, 2*pi)
        X = np.array([
            [r * cos(_theta)],
            [r * sin(_theta)],
            [rand.uniform(-L/2, L/2)],
            ])
        obs.append(X)
    
    return _rotate(theta, phi, zeta) @ obs + center


def _set_field(lx, ly, center, n, theta=0, phi=0, zeta=0):
    """面を表現"""
    
    obs = []
    
    rand = np.random.RandomState(123)
    for i in range(n):
        X = np.array([
            [rand.uniform(-1, 1) * lx/2],
            [rand.uniform(-1, 1) * ly/2],
            [0],
            ])
        obs.append(X)
    
    return _rotate(theta, phi, zeta) @ obs + center


def _set_box(lx, ly, lz, center, n, theta=0, phi=0, zeta=0,):
    
    pass


def set_obstacle(data=None):
    """固定障害物を返す"""
    
    def _choice(name):
        if name == "cylinder":
            return _set_cylinder
        elif name == "sphere":
            return _set_sphere
        elif name == 'field':
            return _set_field
    
    
    if data is None:
        return None
    else:
        obs = []
        for d in data:
            obs.extend(
                _choice(d[0])(*d[1:])
            )
    return obs


data1 = [
    ["cylinder", 0.1, 1.0, np.array([[0.25, -0.4, 1]]).T, 50, 0, pi/2, 0],
    ["cylinder", 0.1, 1.0, np.array([[0.25, -0.4, 1]]).T, 50, pi/2, 0, 0]
]


data2 = [
    ['field', 1, 0.5, np.array([[0.25, -0.5, 1.15]]).T, 100, 0, 0, 0]
]  # 机



def _test(data):

    right = BaxterRobotArmKinematics(isLeft=False)
    os_r = right.get_joint_positions()
    xrs, yrs, zrs = [], [], []
    for o in os_r:
        xrs.append(o[0, 0])
        yrs.append(o[1, 0])
        zrs.append(o[2, 0])

    left = BaxterRobotArmKinematics(isLeft=True)
    os_l = left.get_joint_positions()
    xls, yls, zls = [], [], []
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
    for i, cs in enumerate(right.cpoints_x):
        cs_ = np.concatenate(cs, axis=1)
        xs = cs_[0, :].tolist()
        ys = cs_[1, :].tolist()
        zs = cs_[2, :].tolist()
        ax.scatter(xs, ys, zs, label = "R-" + cs_name[i])
    for i, cs in enumerate(left.cpoints_x):
        cs_ = np.concatenate(cs, axis=1)
        xs = cs_[0, :].tolist()
        ys = cs_[1, :].tolist()
        zs = cs_[2, :].tolist()
        ax.scatter(xs, ys, zs, label = "L-" + cs_name[i])


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



    #ax.legend()
    ax.set_box_aspect((1, 1, 1))



    obs = set_obstacle(data)
    obs = np.concatenate(obs, axis=1)
    ax.scatter(
        obs[0,:], obs[1,:], obs[2,:], marker='.', color = 'k',
    )



    
    plt.show()



if __name__ == '__main__':
    _test(data2)