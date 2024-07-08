from tkinter import NS
import numpy as np
import random

class Bunch:
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

def lam(Nstates):
    
    lambda_1 = (1 + np.sqrt(2 * Nstates - 1)) / 2.0
    lambda_2 = (1 - np.sqrt(2 * Nstates - 1)) / 2.0

    return lambda_1, lambda_2

def initMapping(Nstates, initState = 0):

    delta_1 = np.zeros((Nstates), dtype=float)
    sigma_1 = np.zeros((Nstates), dtype=float)
    for i in range(Nstates):
        delta_1[i] = 2 * (np.random.randint(0,2) - 0.5)
        sigma_1[i] = 2 * (np.random.randint(0,2) - 0.5)

    lambda_1, lambda_2 = lam(Nstates)

    c_1 = np.zeros((Nstates), dtype=complex)
    c_2 = np.zeros((Nstates), dtype=complex)

    factor_1 = np.sqrt( lambda_1**2 / (lambda_1**2 + (Nstates - 1)/2.0) )
    factor_2 = np.sqrt( lambda_2**2 / (lambda_2**2 + (Nstates - 1)/2.0) )
    c_1[0] = factor_1
    c_2[0] = factor_2

    for i in range(1, Nstates):
        c_1[i] = factor_1 * (delta_1[i - 1] + 1.0j * sigma_1[i - 1]) / (2 * lambda_1)
        c_2[i] = factor_2 * (delta_1[i - 1] + 1.0j * sigma_1[i - 1]) / (2 * lambda_2)

    c_1 = np.roll(c_1, initState)
    c_2 = np.roll(c_2, initState)

    return c_1, c_2

def Umap(z, dt, VMat):
        
    Zreal = np.real(z) 
    Zimag = np.imag(z) 

    # Propagate Imaginary first by dt/2
    Zimag -= 0.5 * VMat @ Zreal * dt
    # Propagate Real by full dt
    Zreal += VMat @ Zimag * dt
    # Propagate Imaginary final by dt/2
    Zimag -= 0.5 * VMat @ Zreal * dt

    return  Zreal + 1j*Zimag

def pop(dat):

    c1s, c2s = dat.c1s, dat.c2s
    Nstates = len(c1s)

    lambda_1, lambda_2 = lam(Nstates)

    return lambda_1 * np.outer(c1s, c1s.conjugate()) + lambda_2 * np.outer(c2s, c2s.conjugate())

def Force(dat):

    R  = dat.R
    dH = dat.dHij  
    dH0 = dat.dH0
    NStates = dat.param.NStates

    η = pop(dat)
    η = np.real(η)

    F = np.zeros((len(R)))
    F -= dH0
    for i in range(NStates):
        F -= dH[i,i,:] * η[i,i]
        for j in range(i+1,NStates): # Double counting off-diagonal to save time
            F -= 2.0 * dH[i,j,:] * η[i,j]

    return F

def VelVer(dat) : 

    c1s, c2s = dat.c1s, dat.c2s
    par =  dat.param
    v = dat.P/par.M         
    EStep = int(par.dtN/par.dtE)

    for t in range(EStep):
        c1s = Umap(c1s, par.dtE/2.0, dat.Hij)
        c2s = Umap(c2s, par.dtE/2.0, dat.Hij)
    dat.c1s, dat.c2s = c1s * 1, c2s * 1

    F1    =  Force(dat) # force with {qF(t+dt/2)} * dH(R(t))
    dat.R += v * par.dtN + 0.5 * F1 * par.dtN ** 2 / par.M      
    dat.Hij  = par.Hel(dat.R)
    dat.dHij = par.dHel(dat.R)
    dat.dH0  = par.dHel0(dat.R)
    #-----------------------------
    F2 = Force(dat) # force with {qF(t+dt/2)} * dH(R(t+ dt))
    v += 0.5 * (F1 + F2) * par.dtN / par.M
    dat.P = v * par.M
    dat.Hij = par.Hel(dat.R)

    for t in range(EStep):
        c1s = Umap(c1s, par.dtE/2.0, dat.Hij)
        c2s = Umap(c2s, par.dtE/2.0, dat.Hij)
    dat.c1s, dat.c2s = c1s * 1, c2s * 1

    return dat

def runTraj(parameters):

    try: 
        np.random.seed(parameters.SEED)
    except:
        pass
    #------------------------------------
    NSteps = parameters.NSteps
    NTraj = parameters.NTraj
    NStates = parameters.NStates
    initState = parameters.initState
    nskip = parameters.nskip
    
    pl = 0 if (NSteps % nskip == 0) else 1
    
    rho_ensemble = np.zeros((NStates,NStates,NSteps//nskip + pl), dtype=complex)

    for itraj in range(NTraj):      # loop over trajectories

        dat = Bunch(param =  parameters )

        dat.c1s, dat.c2s = initMapping(NStates, initState)
        dat.R, dat.P = parameters.initR()

        iskip = 0

        dat.Hij  = parameters.Hel(dat.R)
        dat.dHij = parameters.dHel(dat.R)
        dat.dH0  = parameters.dHel0(dat.R)
        
        iskip = 0

        for i in range(NSteps):     
            if (i % nskip == 0):        
                rho_ensemble[:,:,iskip] += pop(dat)
                iskip += 1
            dat = VelVer(dat)

    return rho_ensemble

