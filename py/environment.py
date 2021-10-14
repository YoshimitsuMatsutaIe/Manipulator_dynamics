"""シミュレーション環境"""


import numpy as np

np.random.seed(0)  # 固定

from math import cos, sin, tan, pi
import matplotlib.pyplot as plt


from kinematics import BaxterRobotArmKinematics


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


def _set_point(center):
    """点を置く"""
    return [np.array(center).T]

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
    
    R = _rotate(
        np.deg2rad(theta),
        np.deg2rad(phi),
        np.deg2rad(zeta),
    )
    return R @ obs + np.array(center).T


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


def set_obstacle(obs_param):
    """固定障害物を返す"""
    
    def _choice(name):
        if name == 'point':
            return _set_point
        elif name == "cylinder":
            return _set_cylinder
        elif name == "sphere":
            return _set_sphere
        elif name == 'field':
            return _set_field
    
    
    if obs_param is None:
        return None
    else:
        obs = []
        for d in obs_param:
            obs.extend(
                _choice(d['name'])(**d['data'])
            )
    return obs


data1 =[
    {
        'name' : 'cylinder',
        'data' : {
            'r' : 0.1,
            'L' : 1.0,
            'center' : [[0.25, -0.7, 1]],
            'n' : 100,
            'theta' : 0,
            'phi' : 0,
            'zeta' : 0,
        },
    },
    {
        'name' : 'cylinder',
        'data' : {
            'r' : 0.06,
            'L' : 1.0,
            'center' : [[-0.25, -0.4, 1.25]],
            'n' : 100,
            'theta' : 60,
            'phi' : 30,
            'zeta' : 30,
        },
    },
]




class Goal:

    def __init__(self, **kwargs):
        
        name = kwargs.pop('name')
        
        if name == 'static':
            self.goal = self.static
            self.center = np.array(kwargs.pop('center')).T
        
        elif name == 'tracking_circle':
            self.goal = self.tracking_circle
            self.center = np.array(kwargs.pop('center')).T
            self.r = kwargs.pop('r')
            self.omega = kwargs.pop('omega')
            self.init_theta = kwargs.pop('init_theta')
            self.alpha = kwargs.pop('alpha')
            self.beta = kwargs.pop('beta')
            self.gumma = kwargs.pop('gumma')
        
        return
    
    
    def tracking_circle(self, t):
        R = _rotate(self.alpha, self.beta, self.gumma)
        X = np.array([[
            self.r * np.cos(self.omega*t),
            self.r * np.sin(self.omega * t),
            0,
        ]]).T
        
        return R @ X + self.center


    def static(self, t):
        return self.center


def _test(data):


    goal = np.array([[0.0, -0.5, 1]]).T


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

    ax.scatter(
        goal[0, 0], goal[1, 0], goal[2, 0],
        s = 100, label = 'goal point', marker = '*', color = '#ff7f00', 
        alpha = 1, linewidths = 1.5, edgecolors = 'red')

    obs = set_obstacle(data)
    obs = np.concatenate(obs, axis=1)
    ax.scatter(
        obs[0,:], obs[1,:], obs[2,:], marker='.', color = 'k',
    )

    
    plt.show()



if __name__ == '__main__':
    _test(data1)