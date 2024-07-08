import numpy as np
from numpy import array as A

# ================= global ====================

conv = 27.2114                              # eV / a.u.
fs_to_au = 41.341                           # a.u./fs

# =============================================

class parameters():
                           
    dtN = 0.01 * fs_to_au                   # Nuclear Time Step, a.u.
    T = 500 * fs_to_au                      # Length of Simulation, a.u.
    NSteps = int( T / dtN )                 # 41.341 a.u. / fs
    NTraj = 100
    EStep = 20
    dtE = dtN / EStep
    nskip = 50                              # Plot every {NSkip} nuclear steps

    # MODEL-SPECIFIC ITEMS
    NStates = 2
    ndof = 3
    initState = 1

    # VAR      STATE 1       STATE 2
    Ek    = A([ 3.940/conv,  4.840/conv  ])
    w1    = A([ 0.126/conv,  0.126/conv  ])
    kap1  = A([ 0.037/conv,  -0.254/conv ])
    w6a   = A([ 0.074/conv,  0.074/conv  ])
    kap6a = A([ -0.105/conv, 0.149/conv  ])
    w10a  = A([ 0.118/conv,  0.118/conv  ])

    # auxiliary, dim = len(R)
    omega = A([ w1[0], w6a[0], w10a[0] ])
    kap_1 = A([ kap1[0], kap6a[0], 0.0 ])
    kap_2 = A([ kap1[1], kap6a[1], 0.0 ])
    lam   = A([ 0.0, 0.0, 0.262/conv ])
    M     = A([ 1/w1[0], 1/w6a[0], 1/w10a[0] ]) # 3 Modes = 3 Masses

def Hel(R): # explicitly put as matrix

    # VAR      STATE 1       STATE 2
    Ek      = parameters.Ek
    kap_1   = parameters.kap_1
    kap_2   = parameters.kap_2
    omega   = parameters.omega
    lam     = parameters.lam
    NStates = parameters.NStates
    
    VMat = np.zeros(( NStates, NStates )) # Diabatic Hamiltonian
    
    VMat[0,0] += Ek[0] + 0.5 * np.dot( omega, R**2 ) + np.dot( kap_1, R )
    VMat[1,1] += Ek[1] + 0.5 * np.dot( omega, R**2 ) + np.dot( kap_2, R )

    VMat[0,1] += np.dot( lam, R )
    VMat[1,0] = VMat[0,1] * 1.0

    return VMat

def dHel0(R):

    ndof = parameters.ndof
    dVMat0 = np.zeros(( ndof ))
    w = parameters.omega
    dVMat0 = w * R  # still a matrix

    return dVMat0

def dHel(R):

    kap1    = parameters.kap1
    kap6a   = parameters.kap6a
    lam     = parameters.lam
    NStates = parameters.NStates
    ndof    = parameters.ndof

    dVMat = np.zeros(( NStates, NStates, ndof ))

    dVMat[0,0,0] = kap1[0]  # First Mode
    dVMat[0,0,1] = kap6a[0] # Second Mode

    dVMat[0,1,2] = lam[2]  # Third Mode
    dVMat[1,0,2] = dVMat[0,1,2]

    dVMat[1,1,0] = kap1[1]  # First Mode
    dVMat[1,1,1] = kap6a[1] # Second Mode

    return dVMat

def initR(): # Initial product state of the vibrational ground state <\Psi| = \Product_j \pi^{-1/4} \exp{ -x_j^2 / 2 }

    ndof    = parameters.ndof
#    w       = parameters.omega
    
    R = np.zeros( ndof )
    P = np.zeros( ndof )  

    sigP = 1.0
    sigR = 1.0

    for d in range( ndof ):
        R[d] = np.random.normal() * sigR    #[d]       # np.random.normal() gives exp( -x**2 / 2) distribution
        P[d] = np.random.normal() * sigP    #[d]  
        
    return R, P