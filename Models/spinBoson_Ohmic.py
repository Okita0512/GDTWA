import numpy as np
from numpy import array as A

def model(K):

    #        |     symmetric      |              asymmetric                 | decoupled
    #        |  K0  |  K1  |  K2  |  K3  |  K4  |  K5  |  K6  |  K7  |  K8  |  K9
    ε   = A([  0.00,  0.00,  0.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  0.0  ])         # system energy level
    Δ   = A([  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.0  ])         # coupling
    β   = A([  0.25,  0.25,  0.25,  0.25,  5.00,  5.00,  5.00,  5.00,  10.0,  1.0  ])         # temperature reverse
    ωc  = A([  5.00,  1.00,  0.25,  1.00,  2.00,  2.50,  2.50,  2.50,  2.50,  1.0  ])         # cut-off frequency
    α   = A([  0.02,  0.10,  0.40,  0.40,  0.40,  0.10,  0.20,  0.40,  0.20,  0.0  ])         # Kondo parameter
    N   = A([  50,    50,    100,   100,   50,    50,    50,    100,   50,    1  ])           # number of bath oscillators

    return ε[K], Δ[K], β[K], ωc[K], α[K], N[K]  

# Symmetric: Jian Liu's 2018 paper Model 3 - 5
# Asymmetric: Jian Liu's 2018 paper Model 6 - 11

def bathParam(ωc, α, ndof):     # bath descretization

#    ωm = 4.0            # cut-off frequency that is sufficiently larger than ωc
#    ω0 = ωc * ( 1 - np.exp(- ωm) ) / ndof       # should be ωm/ωc ?

    c = np.zeros(ndof)
    ω = np.zeros(ndof)
    for d in range(ndof):
        """
        2018 Jian Liu, 2005 H. Wang and S. Thoss
        """
        ω[d] =  - ωc * np.log(1 - (d + 1)/(ndof + 1))
        c[d] =  np.sqrt(α * ωc / ((ndof + 1))) * ω[d]

        """
        2004 Coker
        """
#        c[d] =  np.sqrt(α * ω0) * ω[d] 

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 1D spin boson test !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#    c = 1.0
#    ω = 1.0
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    return c, ω


class parameters():
   NSteps = 10     # total time = NSteps * dtN (a.u.)
   NTraj = 1      # number of trajectories, 10**6 for publication
   dtN = 0.01       # Nuclei step, low temperature use 0.002. Here for all these models 0.01 is enough.
   dtE = dtN/10     # Electronic step
   NStates = 2      # spin-boson
   M = 1            # Mass, atomic unit 1 for electron
   K = 5            # model K
   initState = 0    #
   nskip = 10
   ε, Δ, β, ωc, α, ndof = model(K)
#   ε, Δ, β, ωc, α, ndof =        
   c, ω  = bathParam(ωc, α, ndof)

# for NRPMD only
   nb = 8       # kill this. make nb to be included in main program


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
#    c = parameters.c
    ω = parameters.ω

    dH0 = ω**2 * R 

    return dH0


def dHel(R):
    c = parameters.c
#    ω = parameters.ω
    
    dHij = np.zeros(( 2,2,len(R)))
    dHij[0,0,:] = c   
    dHij[1,1,:] = - c         
    return dHij         

def initR():            # initial state of the bath
    R0 = 0.0            
    P0 = 0.0
    β  = parameters.β
    ω = parameters.ω
    ndof = parameters.ndof

# obtained from Wigner transform of the bath probability density

# choice 1: bath is initially in thermal equilibrium 标注下文献来源

    sigP = np.sqrt( ω / ( 2 * np.tanh( 0.5*β*ω ) ) )
    sigR = sigP / ω

# choice 2: the canonical density for the nuclei part 标注下文献来源

#    sigP = np.sqrt( ω / ( 2 * np.tanh( 0.5*β*ω ) ) )
#    sigR = (sigP + c / ω**2 ) / ω

    R = np.zeros(( ndof ))
    P = np.zeros(( ndof ))
    for d in range(ndof):
        R[d] = np.random.normal()*sigR[d]
        P[d] = np.random.normal()*sigP[d]  
    return R, P

# for 1D only
#    R = np.random.normal()*sigR
#    P = np.random.normal()*sigP
#    return R, P
