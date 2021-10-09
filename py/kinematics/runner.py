"""メイン"""

import yaml

import rmp_simulation
#import param





def run(params):
    
    with open(params) as file:
        config = yaml.safe_load(file.read())
    
    sim_param = config['sim_param']
    rmp_param = config['rmp_param']
    
    simulator = rmp_simulation.Simulator(**sim_param)
    simulator.set_controller(rmp_param)
    
    
    return



if __name__ == '__main__':
    #run(*param.setting_1())
    
    run('./py/kinematics/sean_1.yaml')
