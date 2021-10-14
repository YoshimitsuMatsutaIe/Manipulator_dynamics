"""メイン"""

import yaml


import rmp_simulation




def run(params):
    """シミュレーションを実行
    
    params : yamlでシミュレーション条件を教えて
    """
    
    with open(params, encoding='UTF-8') as file:
        config = yaml.safe_load(file.read())
    
    sim_param = config['sim_param']
    rmp_param = config['rmp_param']
    env_param = config['env_param']
    
    
    simulator = rmp_simulation.Simulator(**sim_param)
    simulator.set_controller(rmp_param)
    simulator.set_environment(env_param)
    simulator.run_simulation()
    simulator.plot_animation_2()
    
    return


if __name__ == '__main__':
    
    #run('./config/scean_1.yaml')
    
    #run('./config/tracking.yaml')
    
    run('./config/use_RMPfromGDS_test.yaml')
    
    #run('.//config/test.yaml')
