import numpy as np
from math import pi

from baxter_utils_3 import BaxterKinematics3
from old import BaxterKinematics



q = np.array([[0, -31, 0, 43, 0, 72, 0]]).T * pi/180 

hoge = BaxterKinematics3()
J_all = hoge.jacobi_all(q)
print(len(J_all))

hogehoge = BaxterKinematics()
J_all_2 = hogehoge.Jo_global_l
print(len(J_all_2))

# d = J_all[2] - J_all_2[0]

# print(d)



print(J_all[4])
print(J_all_2[0])