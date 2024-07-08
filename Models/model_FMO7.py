import numpy as np
from numpy import random as rn
from numba import jit

au = 27.2113961 # 27.2114
ps = 41341.374575751
fs = 41.34136 # 41.341374575751
cm = 1/4.556335e-6 # 219474.63
kB = 3.166829e-6

NSteps = 16000 # number of nuclear steps in the simulation
NStepsPrint = 500 # number of nuclear steps that are stored/printed/outputed
NSkip = int(NSteps/NStepsPrint) # used to enforce NStepsPrint
totalTime = 4000*fs # total amount of simulation time
dtN = totalTime/NSteps # nuclear timestep 
dtE = dtN # /1 # electronic timestep (should be smaller than nuclear)
NTraj = 100
NStates = 7
NBath = 60
NR = NBath*NStates
M = np.ones(NR)
initBasis = 1 # 0 = adiabatic, 1 = exciton (no nuclei), 2 = diabatic
initState = 0 # Site 1
outputBasis = 1 # 0 = adiabatic, 1 = exciton (no nuclei), 2 = diabatic

# MASH-specific parameters

sample_type = 'focused' # 'focused' or 'gaussian', probably not important
hop_type = 'express' # 'express' is faster, 'full' uses max_hop and num_bisect
max_hop = 30
num_bisect = 10

λ = 35/cm
τc = 50*fs
ωc = 106.14/cm # 1/τc
ωk = ωc * np.tan(0.5*np.pi*(np.arange(NBath)+0.5)/NBath)
#ωk = ωc * np.tan( np.pi * (1 - np.arange(1,NBath+1)/(NBath + 1)) / 2)
κ = ωk * np.sqrt(2*λ/NBath)
#κ = ωk * np.sqrt(2*λ/(NBath+1))
β = 1 / (300 * kB) # 1053

@jit(nopython=True)
def Hel(R):
    H = np.array([[200,     -87.7,  5.5,    -5.9,   6.7,    -13.7,  -9.9    ],
                  [-87.7,   320,    30.8,   8.2,    0.7,    11.8,   4.3     ],
                  [5.5,     30.8,   0,      -53.5,  -2.2,   -9.6,   6.0     ],
                  [-5.9,    8.2,    -53.5,  110,    -70.7,  -17.0,  -63.3   ],
                  [6.7,     0.7,    -2.2,   -70.7,  270,    81.1,   -1.3    ],
                  [-13.7,   11.8,   -9.6,   -17.0,  81.1,   420,    39.7    ],
                  [-9.9,    4.3,    6.0,    -63.3,  -1.3,   39.7,   230     ]], dtype = np.complex_)/cm
    for j in range(NStates):
        H[j, j] += np.sum(κ * R[j*NBath:(j+1)*NBath])
    return H

@jit(nopython=True)
def H_const():
    H = np.array([[200,     -87.7,  5.5,    -5.9,   6.7,    -13.7,  -9.9    ],
                  [-87.7,   320,    30.8,   8.2,    0.7,    11.8,   4.3     ],
                  [5.5,     30.8,   0,      -53.5,  -2.2,   -9.6,   6.0     ],
                  [-5.9,    8.2,    -53.5,  110,    -70.7,  -17.0,  -63.3   ],
                  [6.7,     0.7,    -2.2,   -70.7,  270,    81.1,   -1.3    ],
                  [-13.7,   11.8,   -9.6,   -17.0,  81.1,   420,    39.7    ],
                  [-9.9,    4.3,    6.0,    -63.3,  -1.3,   39.7,   230     ]], dtype = np.complex_)/cm
    return H

@jit(nopython=True)
def dHel(R):
    dH = np.zeros((NStates,NStates,NR), dtype = np.complex_)
    for j in range(NStates):
        dH[j,j,j*NBath:(j+1)*NBath] = κ
    return dH

@jit(nopython=True)
def dHel0(R):
    dH0 = np.zeros(NR) # state independent, so only need to do NR derivatives once
    for j in range(NStates):
        dH0[j * NBath : (j + 1) * NBath] = ωk**2 * R[j * NBath : (j + 1) * NBath]
    return dH0

@jit(nopython=True)
def initR():
    R0 = 0.0 # average initial nuclear position (could be array specific to each nuclear DOF)
    P0 = 0.0 # average initial nuclear momentum (could be array specific to each nuclear DOF)
    #σP = np.sqrt(ωk/(2.0 * np.tanh(0.5 * β * ωk))) # standard dev of nuclear momentum (mass = 1)
    #σR = σP/(ωk) # standard dev of nuclear position (mass = 1)
    σP = np.sqrt(1/β)*np.ones(NBath)
    σR = 1/np.sqrt(β*ωk**2)
    R = np.zeros((NR))
    P = np.zeros((NR))
    for Ri in range(NStates*NBath):
        imode = (Ri%NBath) # current bath mode
        R[Ri] = rn.normal() * σR[imode] + R0 # gaussian random variable centered around R0
        P[Ri] = rn.normal() * σP[imode] + P0 # gaussian random variable centered around P0
    return R, P