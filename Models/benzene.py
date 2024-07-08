import numpy as np
from numpy import array as A

# ================= global ====================

conv = 27.211397                              # eV / a.u.
fs_to_au = 41.341                           # a.u./fs

class parameters():
                           
    dtN = 0.01 * fs_to_au                   # Nuclear Time Step, a.u.
    T = 1 * fs_to_au                      # Length of Simulation, a.u.
    NSteps = int( T / dtN )                 # 41.341 a.u. / fs
    NTraj = 1
    EStep = 20
    dtE = dtN / EStep
    nskip = 50                              # Plot every {NSkip} nuclear steps

    # MODEL-SPECIFIC ITEMS
    NStates = 3
    ndof = 5
    initState = 2

    # VAR      STATE 1       STATE 2          STATE 3
    Ek    = A([ 9.75/conv,    11.84/conv,     12.44/conv  ])

    # auxiliary, dim = len(R)
    kap_1  = A([ -0.042/conv, -0.246/conv, -0.125/conv, 0.0, 0.0 ])
    kap_2  = A([ -0.042/conv, 0.242/conv,  0.1/conv,    0.0, 0.0 ])
    kap_3  = A([ -0.301/conv, 0.0,         0.0,         0.0, 0.0 ])
    omega    = A([ 0.123/conv, 0.198/conv, 0.075/conv, 0.088/conv, 0.12/conv ])
    lam_12   = A([ 0.0, 0.0, 0.0, 0.164/conv, 0.0 ])
    lam_23   = A([ 0.0, 0.0, 0.0, 0.0, 0.154/conv ])    # no coupling between 1 and 3

    # auxiliary 
    M = A([ 1/omega[0], 1/omega[1], 1/omega[2], 1/omega[3], 1/omega[4] ]) # 5 Modes = 5 Masses

def Hel(R):

    Ek      = parameters.Ek
    kap_1   = parameters.kap_1
    kap_2   = parameters.kap_2
    kap_3   = parameters.kap_3
    omega   = parameters.omega
    lam_12     = parameters.lam_12
    lam_23     = parameters.lam_23
    NStates = parameters.NStates
    
    VMat = np.zeros(( NStates, NStates )) # Diabatic Hamiltonian
    
    VMat[0,0] += Ek[0] + 0.5 * np.dot( omega, R**2 ) + np.dot( kap_1, R )
    VMat[1,1] += Ek[1] + 0.5 * np.dot( omega, R**2 ) + np.dot( kap_2, R )
    VMat[2,2] += Ek[2] + 0.5 * np.dot( omega, R**2 ) + np.dot( kap_3, R )

    VMat[0,1] += np.dot( lam_12, R )
    VMat[1,0] = VMat[0,1] * 1.0

    VMat[1,2] += np.dot( lam_23, R )
    VMat[2,1] = VMat[1,2] * 1.0

    return VMat

def dHel0(R):

    ndof = parameters.ndof
    dVMat0 = np.zeros(( ndof ))
    w = parameters.omega
    dVMat0 = w * R

    return dVMat0

def dHel(R):

    kap_1   = parameters.kap_1
    kap_2   = parameters.kap_2
    kap_3   = parameters.kap_3
    lam_12  = parameters.lam_12
    lam_23  = parameters.lam_23
    NStates = parameters.NStates
    ndof    = parameters.ndof

    dVMat = np.zeros(( NStates, NStates, ndof ))

    # diagonal coupling
    dVMat[0,0,0] = kap_1[0]  # First Mode
    dVMat[0,0,1] = kap_1[1]  # Second Mode
    dVMat[0,0,2] = kap_1[2]
    dVMat[1,1,0] = kap_2[0]  # First Mode
    dVMat[1,1,1] = kap_2[1]  # Second Mode
    dVMat[1,1,2] = kap_2[2]
    dVMat[2,2,0] = kap_3[0]

    # off-diagonal coupling
    dVMat[0,1,3] = lam_12[3]
    dVMat[1,0,3] = dVMat[0,1,3]
    dVMat[1,2,4] = lam_23[4]
    dVMat[2,1,4] = dVMat[1,2,4]

    return dVMat

def initR(): # Initial product state of the vibrational ground state <\Psi| = \Product_j \pi^{-1/4} \exp{ -x_j^2 / 2 }

    ndof    = parameters.ndof
    
    R = np.zeros( ndof )
    P = np.zeros( ndof )  

    sigP = 1.0
    sigR = 1.0

    for d in range( ndof ):
        R[d] = np.random.normal() * sigR    #[d]       # np.random.normal() gives exp( -x**2 / 2) distribution
        P[d] = np.random.normal() * sigP    #[d]  
        
    return R, P

# print(Hel(np.array([1,1,1,1,1])))