"""メイン"""

import rmp_simulation
import param





def run(sim_param, rmp_param, env_param):
    
    simulator = rmp_simulation.Simulator(**sim_param)
    
    simulator.set_controller(rmp_param)
    
    return



if __name__ == '__main__':
    run(*param.setting_1())
