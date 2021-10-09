import numpy as np
from math import cos, sin, tan, pi

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
    rotate = np.array([
        [1, 0, 0],
        [0, cos(theta), -sin(theta)],
        [0, sin(theta), cos(theta)],
    ]) @ np.array([
        [cos(phi), 0, sin(phi)],
        [0, 1, 0],
        [-sin(phi), 0, cos(phi)],
    ]) @ np.array([
        [cos(zeta), -sin(zeta), 0],
        [sin(zeta), cos(zeta), 0],
        [0, 0, 1],
    ])
    
    rand = np.random.RandomState(123)
    for i in range(n):
        _theta = rand.uniform(0, 2*pi)
        X = np.array([
            [r * cos(_theta)],
            [r * sin(_theta)],
            [rand.uniform(-L/2, L/2)],
            ])
        obs.append(rotate @ X + center)
    
    return obs




def set_obstacle(data=None):
    """固定障害物を返す"""
    
    def _choice(name):
        if name == "cylinder":
            return _set_cylinder
        elif name == "sphere":
            return _set_sphere
    
    
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
    ["cylinder", 0.1, 1.0, np.array([[0.25, -0.4, 1]]).T, 5, 0, 0, 0],
    ["cylinder", 0.1, 1.0, np.array([[0.25, -0.4, 1]]).T, 5, 0, pi/2, 0]
]






if __name__ == '__main__':
    
    o = set_obstacle(data1)
    print(len(o))