import numpy as np

def centered_derivative( time, com) :
    """Centered finite difference of com wrt time

    the first point is computed with forward and the last with backward differences
    """
    ext_time = time
    ext_time = np.append( ext_time, 2*time[-1]-time[-2])
    ext_time = np.append( 2*time[0]-time[1], ext_time)
    ext_com = com
    ext_com = np.append( ext_com, 2*com[-1]-com[-2])
    ext_com = np.append( 2*com[0]-com[1], ext_com)
    centered = [-1,0,1]
    return np.convolve( centered, ext_com, 'valid')/ np.convolve(
            centered, ext_time, 'valid')

def generate_default( kappa = 0.000457, amplitude = 1, sigma = 10):
    """Default real simulation parameters"""
    return {
        "grid" :
        {
            "n"  : 4,
            "Nx" : 512,
            "Ny" : 512,
            "x" : [0, 500],
            "y" : [0, 500]
        },
        "output":
        {
            "type" : "netcdf",
            "tend" : 1500,
            "maxout" : 10,
            "itstp": 30,
            "n"  : 4,
            "Nx" : 256,
            "Ny" : 256
        },
        "elliptic":
        {
            "stages" : 4,
            "eps_pol" : [1e-6, 0.5, 0.5,0.5],
            "eps_gamma" : [1e-8, 1, 1,1],
            "direction" : "forward"
        },
        "timestepper":
        {
            "type" : "adaptive",
            "tableau" : "Bogacki-Shampine-4-2-3",
            "rtol" : 1e-6,
            "atol" : 1e-7,
            "reject-limit" : 2
        },
        "curvature" : kappa,
        "potential":
        {
            "epsilon_D" : 0,
            "bc" : ["DIR", "PER"]
        },
        "species":
        [
            {
                "name" : "e",
                "mu" : 0,
                "tau" : -1,
                "a" : -1,
                "nu_perp" : 1e-8, # 1e-5
                "bc" : ["DIR", "PER"],
                "init":
                {
                    "type" : "zero_potential"
                }
            },
            {
                "name" : "i",
                "mu" : 1,
                "tau" : 0.0,
                "a" : 0.9,
                "nu_perp" : 1e-8, # 1e-5
                "bc" : ["DIR", "PER"],
                "init":
                {
                    "type" : "blob",
                    "amplitude" : amplitude,
                    "posX" : 0.5,
                    "posY" : 0.5,
                    "sigma" : sigma,
                    "flr" : "gamma_inv"
                },
            },
            {
                "name" : "j",
                "mu" : 2,
                "tau" : 0.0,
                "a" : 0.1,
                "nu_perp" : 1e-8, # 1e-5
                "bc" : ["DIR", "PER"],
                "init":
                {
                    "type" : "blob",
                    "amplitude" : amplitude,
                    "posX" : 0.5,
                    "posY" : 0.5,
                    "sigma" : sigma,
                    "flr" : "gamma_inv"
                }
            }
        ]
    }

