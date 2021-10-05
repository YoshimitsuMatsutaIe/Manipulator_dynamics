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
    



def main():
    arm = BaxterRobotArmKinematics(isLeft=True)
    




if __name__ == "__main__":
    print("実行中...")
    start = time.time()
    
    main()
    
    print("実行終了")
    print("実行時間 = ", time.time() - start)