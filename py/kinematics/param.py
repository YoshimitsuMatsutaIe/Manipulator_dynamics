"""パラメータ"""


def setting_1():
    """パターン1
    静止目標（障害物なし）
    
    """
    sim_param = {
        'isleft' : True,
        'TIME_SPAN' : 10,
        'TIME_INTERVAL' : 0.05,
    }
    rmp_param = [
        {
            'name' : 'original',
            
        }
    ]
    
    env_param = {
        
    }
    
    return sim_param, rmp_param, env_param