import simplesimdb as simplesim
import numpy as np
import json

kappa = 0.000457
amplitude = 1
sigma = 10

inputfile={
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

m = simplesim.Manager( directory="vmax", executable="./submit_job.sh", filetype="nc")

tau = (0,0)
# This is for Figure 2 in the paper
for aa in (1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 1e-1, 2e-1, 5e-1):
    for mumu in ( 2, 2.5, 5, 10, 20):
#for aa in [5e-1]:
#    for mumu in [ 20]:
        a = ( 1-aa, aa)
        mu = ( 1, mumu)
        vmax = 0.77*np.sqrt( amplitude*sigma*kappa*(1+a[0]*tau[0]+a[1]*tau[1])/(a[0]*mu[0]+a[1]*mu[1]) )
        inputfile["species"][1]["mu"] = mu[0]
        inputfile["species"][2]["mu"] = mu[1]
        inputfile["species"][1]["a"] = a[0]
        inputfile["species"][2]["a"] = a[1]
        inputfile["output"]["tend"] = 2*sigma/vmax/0.2

        print ( a[1], mu[1], vmax, sigma/vmax/0.2);
        m.create( inputfile, error="display")


############################################
print( "Figure 4")
m = simplesim.Manager( directory="caseA", executable="./submit_job.sh", filetype="nc")

# Figure 4
aa = 0.01
mumu = 2
a = ( 1-aa, aa)
mu = ( 1, mumu)
tau = ( 2,2)

vmax = 0.77*np.sqrt( amplitude*sigma*kappa*(1+a[0]*tau[0]+a[1]*tau[1])/(a[0]*mu[0]+a[1]*mu[1]) )
inputfile["species"][1]["mu"] = mu[0]
inputfile["species"][2]["mu"] = mu[1]
inputfile["species"][1]["a"] = a[0]
inputfile["species"][2]["a"] = a[1]
inputfile["species"][1]["tau"] = tau[0]
inputfile["species"][2]["tau"] = tau[1]
inputfile["output"]["tend"] = 2*sigma/vmax/0.2

print ( a[1], mu[1], vmax, sigma/vmax/0.2);
m.create( inputfile, error="display")
############################################
print( "Simulation case A")

aa = 0.001
mumu = 2
a = ( 1-aa, aa)
mu = ( 1, mumu)
tau = ( 1,1)

vmax = 0.77*np.sqrt( amplitude*sigma*kappa*(1+a[0]*tau[0]+a[1]*tau[1])/(a[0]*mu[0]+a[1]*mu[1]) )
inputfile["species"][1]["mu"] = mu[0]
inputfile["species"][2]["mu"] = mu[1]
inputfile["species"][1]["a"] = a[0]
inputfile["species"][2]["a"] = a[1]
inputfile["species"][1]["tau"] = tau[0]
inputfile["species"][2]["tau"] = tau[1]
inputfile["output"]["tend"] = 2*sigma/vmax/0.2
inputfile["species"][1]["init"] ={ 
    "type" : "blob",
    "amplitude" : 1.0,
    "posX" : 0.4,
    "posY" : 0.5,
    "sigma" : 10.,
    "flr" : "gamma_inv"
}
inputfile["species"][2]["init"] ={
    "type" : "wall",
    "amplitude" : 5.0,
    "posX" : 0.5,
    "sigma" : 10.,
    "flr" : "gamma_inv"
}
print ( a[1], mu[1], vmax, sigma/vmax/0.2);
m.create( inputfile, error="display")

############################################
print( "Simulation case B")
m = simplesim.Manager( directory="caseB", executable="./submit_job.sh", filetype="nc")
aa = 0.1
mumu = 5 # macht sonst kein Sinn
a = ( 1-aa, aa)
mu = ( 1, mumu)
tau = ( 0,0)

vmax = 0.77*np.sqrt( amplitude*sigma*kappa*(1+a[0]*tau[0]+a[1]*tau[1])/(a[0]*mu[0]+a[1]*mu[1]) )
inputfile["species"][1]["mu"] = mu[0]
inputfile["species"][2]["mu"] = mu[1]
inputfile["species"][1]["a"] = a[0]
inputfile["species"][2]["a"] = a[1]
inputfile["species"][1]["tau"] = tau[0]
inputfile["species"][2]["tau"] = tau[1]
inputfile["output"]["tend"] = 2*sigma/vmax/0.2
inputfile["species"][1]["init"] ={ 
    "type" : "blob",
    "amplitude" : 1.0,
    "posX" : 0.3,
    "posY" : 0.5,
    "sigma" : 10.,
    "flr" : "gamma_inv"
}
inputfile["species"][2]["init"] ={
    "type" : "wall",
    "amplitude" : 5.0,
    "posX" : 0.45,
    "sigma" : 10.,
    "flr" : "gamma_inv"
}
print ( a[1], mu[1], vmax, sigma/vmax/0.2);
m.create( inputfile, error="display")

############################################
print( "Simulation case C")
m = simplesim.Manager( directory="caseC", executable="./submit_job.sh", filetype="nc")
tau = (0,0)
# This is for Figure 6 in the paper
for aa in (0, 0.03, 0.075, 0.1, 0.15, 0.2):
    for mumu in ( 1, 3, 5, 7):
        a = ( 1-aa, aa)
        mu = ( 1, mumu)
        vmax = 0.77*np.sqrt( amplitude*sigma*kappa*(1+a[0]*tau[0]+a[1]*tau[1])/(a[0]*mu[0]+a[1]*mu[1]) )
        inputfile["species"][1]["mu"] = mu[0]
        inputfile["species"][2]["mu"] = mu[1]
        inputfile["species"][1]["a"] = a[0]
        inputfile["species"][2]["a"] = a[1]
        inputfile["species"][1]["tau"] = tau[0]
        inputfile["species"][2]["tau"] = tau[1]
        inputfile["output"]["tend"] = 2*sigma/vmax/0.2
        inputfile["species"][1]["init"] ={ 
            "type" : "blob",
            "amplitude" : 1.0,
            "posX" : 0.3,
            "posY" : 0.5,
            "sigma" : 10.,
            "flr" : "gamma_inv"
        }
        inputfile["species"][2]["init"] ={
            "type" : "wall",
            "amplitude" : 5.0,
            "posX" : 0.3,
            "sigma" : 10.,
            "flr" : "gamma_inv"
        }
        print ( a[1], mu[1], vmax, sigma/vmax/0.2);
        m.create( inputfile, error="display")

