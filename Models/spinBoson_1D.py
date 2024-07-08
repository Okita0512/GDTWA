import numpy as np
from numpy import array as A

def model(K):

    ε   = A([  0.0  ])         # system energy level
    Δ   = A([  0.5  ])         # coupling
    β   = A([  16.0  ])        # temperature reverse
    ωc  = A([  1.0  ])         # cut-off frequency
    α   = A([  0.0  ])         # Kondo parameter
    N   = A([  1  ])           # number of bath oscillators

    return ε[K], Δ[K], β[K], ωc[K], α[K], N[K]  

def bathParam():

#    c = A([ np.sqrt(2.0) / 10 ])       # Model A
    c = A([np.sqrt(2.0) / 2])        # Model B
#    c = A([np.sqrt(2.0)])            # Model C
    ω = A([1.0])

    return c, ω

class parameters():
   NSteps = 1    # total time = NSteps * dtN (a.u.)
   NTraj = 5000       # number of trajectories, 10**6 for publication
   dtN = 0.01       # Nuclei step, low temperature use 0.002. Here for all these models 0.01 is enough.
   dtE = dtN/20     # Electronic step
   NStates = 2      # spin-boson
   M = 1            # Mass, atomic unit 1 for electron
   K = 0            # model K
   initState = 0    #
   nskip = 10
   ε, Δ, β, ωc, α, ndof = model(K)
#   ε, Δ, β, ωc, α, ndof =        
   c, ω  = bathParam()

def Hel(R):

    c = parameters.c
    Δ = parameters.Δ 
    ε = parameters.ε

    Vij = np.zeros((2,2))

    Vij[0,0] =   np.sum(c * R) + ε
    Vij[1,1] = - np.sum(c * R) - ε
    Vij[0,1], Vij[1,0] = Δ, Δ 

    return Vij


def dHel0(R):

    ω = parameters.ω
    dH0 = ω**2 * R 

    return dH0

def dHel(R):
    c = parameters.c
#    ω = parameters.ω
    
    dHij = np.zeros((2,2,1))
    dHij[0,0,:] = c   
    dHij[1,1,:] = - c         

    return dHij         

# =================================================================================================
# =================================================================================================

def initR():        # useless, but needed to avoid bugs
    
    R0 = 0.0            
    P0 = 0.0
    β  = parameters.β
    ω = parameters.ω
    ndof = parameters.ndof

# obtained from Wigner transform of the bath probability density

# choice 1: bath is initially in thermal equilibrium
    sigP = np.sqrt( ω / ( 2 * np.tanh( 0.5*β*ω ) ) )
    sigR = sigP / ω

# choice 2: the canonical density for the nuclei part
#    sigP = np.sqrt( ω / ( 2 * np.tanh( 0.5*β*ω ) ) )
#    sigR = (sigP + c / ω**2 ) / ω

    R = np.zeros(( ndof ))
    P = np.zeros(( ndof ))
#    for d in range(ndof):
#        R[d] = np.random.normal()*sigR[d]
#        P[d] = np.random.normal()*sigP[d]  
#    return R, P

# for 1D only
    R = np.random.normal()*sigR
    P = np.random.normal()*sigP
    return R, P
