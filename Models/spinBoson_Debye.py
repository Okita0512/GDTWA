import numpy as np
from numpy import array as A, pi

def model(K):

    #        |             symmetric              |     asymmetric     | decoupled
    #        |   K0   |  K1  |  K2  |  K3  |  K4  |  K5  |  K6  |  K7  |  K8
    ε   = A([  0.000,   0.00,  0.00,  0.00,  0.00,  1.00,  1.00,  1.00,  1.00  ])         # system energy level
    Δ   = A([  1.000,   1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00,  1.00  ])         # coupling
    β   = A([  0.500,   0.50,  5.00,  0.50,  50.0,  0.50,  0.50,  50.0,  1.00  ])         # temperature reverse
    ωc  = A([  0.250,   0.25,  0.25,  5.00,  5.00,  0.25,  5.00,  5.00,  1.00  ])         # cut-off frequency
    λ   = A([  0.025,   0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.00 ])          # coupling strength
    N   = A([  100,     100,   100,   100,   100,   100,   100,   100,   1  ])            # number of bath oscillators

    return ε[K], Δ[K], β[K], ωc[K], λ[K], N[K] 

# Symmetric: Jian Liu's 2018 paper Model 12 - 16
# Asymmetric: Jian Liu's 2018 paper Model 17 - 19

def bathParam(λ, ωc, ndof):

    c = np.zeros(( ndof ))
    ω = np.zeros(( ndof ))
    for d in range(ndof):
        ω[d] =  ωc * np.tan( pi * (1 - (d + 1)/(ndof + 1)) / 2)

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#        c[d] =  np.sqrt(2 * λ) * ω[d]
        c[d] =  np.sqrt(2 * λ / (ndof + 1)) * ω[d]

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    return c, ω

#   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#         propagation time = 30.0 a.u.     ( 1 fs = 242 a.u. )
#   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#     Temp. and coup.    NSteps      dtN         nskip    
#     Low and strong     3000        0.01        10
#     High and weak      15000       0.002       50

#   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

class parameters():
   NSteps = 5000     # 3000 total time = NSteps * dtN (a.u.)
   NTraj = 200      # number of trajectories, 10**6 for publication
   dtN = 0.002       # Nuclei step. Here for all these models use dt = 0.002
   dtE = dtN/20     # Electronic step
   NStates = 2      
   M = 1            
   K = 1            # model K
   initState = 0    
   nskip = 10       # 10 for dtN
   ε, Δ, β, ωc, λ, ndof = model(K)
#   ε, Δ, β, ωc, λ, ndof =          
   c, ω  = bathParam(λ, ωc, ndof)

def Hel(R):
    c = parameters.c
    Δ = parameters.Δ 
    ε = parameters.ε

    Vij = np.zeros((2,2))

    Vij[0,0] =   np.sum(c * R) + ε
    Vij[1,1] = - np.sum(c * R) - ε

    Vij[0,1], Vij[1,0] = Δ, Δ 
    return Vij


def dHel0(R):           #system - bath
    c = parameters.c
    ω = parameters.ω

    dH0 = ω**2 * R 

    return dH0


def dHel(R):            # system - bath
    c = parameters.c
    ω = parameters.ω
    
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

# choice 1: bath is initially in thermal equilibrium

    sigP = np.sqrt( ω / ( 2 * np.tanh( 0.5*β*ω ) ) )
    sigR = sigP / ω         

# choice 2: the canonical density for the nuclei part

#    sigP = np.sqrt( ω / ( 2 * np.tanh( 0.5*β*ω ) ) )
#    sigR = (sigP + c / ω**2 ) / ω

    R = np.zeros(( ndof ))
    P = np.zeros(( ndof ))
    for d in range(ndof):
        R[d] = np.random.normal()*sigR[d]      
        P[d] = np.random.normal()*sigP[d]  
    return R, P

#------ only required for NRPMD----------

#def initHel0(R):
#    M = parameters.M
#    ω = parameters.ω
#    R0 = 0.0
#    
#    return  np.sum(0.5 *M* ω**2 * (R-R0)**2.0)


#------ only required for HEOM ----------

def Hsys():
    Δ = parameters.Δ 
    ε = parameters.ε
    KK = parameters.KK        # Truncation level

    Vij = np.zeros((2,2))

    Vij[0,0] = ε
    Vij[1,1] = - ε

    Vij[0,1], Vij[1,0] = Δ, Δ 
    return Vij