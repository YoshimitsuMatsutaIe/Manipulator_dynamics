import numpy as np
from math import pi

from baxter_utils_3 import BaxterKinematics3
from new import BaxterRobotArmKinematics

import time

q = np.array([[0, -31, 0, 43, 0, 72, 0]]).T * pi/180 

start = time.time()
hoge = BaxterKinematics3()
J_all = hoge.jacobi_all(q)
print("old = ", time.time() - start)
print(len(J_all))

start = time.time()
hogehoge = BaxterRobotArmKinematics(isLeft=True)
J_all_2 = hogehoge.Jo_global
print("new = ", time.time() - start)
print(len(J_all_2))


# d = J_all[-1] - J_all_2[-1]
# print(d)


# print("q0")
# print(J_all[4])
# print(J_all_2[0])

# print("qGL")
# print(J_all[-1])
# print(J_all_2[-1])


for i in range(8):
    print("i = ", i)
    print(J_all[i+3].astype(np.float16))
    print(J_all_2[i].astype(np.float16))

for i in range(8):
    print(i, " = ", np.linalg.norm(J_all[i+3] - J_all_2[i], 2))