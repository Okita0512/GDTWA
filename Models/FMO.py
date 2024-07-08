from tkinter import NS
import numpy as np

# ================= global ====================

conv = 27.211397                            # 1 a.u. = 27.211397 eV
fs_to_au = 41.341                           # 1 fs = 41.341 a.u.
cm_to_au = 4.556335e-06                     # 1 cm^-1 = 4.556335e-06 a.u.
au_to_K = 3.1577464e+05                     # 1 au = 3.1577464e+05 K

def bathParam(λ, ωc, ndof):

    c = np.zeros(( ndof ))
    ω = np.zeros(( ndof ))
    for d in range(ndof):
        ω[d] =  ωc * np.tan( np.pi * (1 - (d + 1)/(ndof + 1)) / 2)
        c[d] =  np.sqrt(2 * λ / (ndof + 1)) * ω[d]

    return c, ω

class parameters():

    dtN = fs_to_au                         # Nuclear Time Step, a.u.
    NSteps = 1000
    NTraj = 100
    EStep = 20
    dtE = dtN / EStep
    nskip = 10                             # Plot every {NSkip} nuclear steps

    # MODEL-SPECIFIC ITEMS
    NStates = 7
    F = 60                                 # number of nuclear DOFs coupling to each state
    ndof = F * NStates                     # total number of nuclear DOFs
    initState = 0
    M = 1
    temperature = 77 / au_to_K
    β = 1 / temperature
    λ = 35 * cm_to_au
    ωc = 106.14 * cm_to_au

    c = np.zeros((NStates, F), dtype=float)
    ω = np.zeros((NStates, F), dtype=float)
    for i in range(NStates):
        c[i, :], ω[i, :]  = bathParam(λ, ωc, F)

def Hel(R):

    c = parameters.c
    F = parameters.F
    NStates = parameters.NStates
    VMat = np.zeros(( NStates, NStates )) # Diabatic Hamiltonian
    
    VMat[0,:] = np.array([12410.0, -87.7, 5.5, -5.9, 6.7, -13.7, -9.9])
    VMat[1,:] = np.array([-87.7, 12530.0, 30.8, 8.2, 0.7, 11.8, 4.3])
    VMat[2,:] = np.array([5.5, 30.8, 12210.0, -53.5, -2.2, -9.6, 6.0])
    VMat[3,:] = np.array([-5.9, 8.2, -53.5, 12320.0, -70.7, -17.0, -63.3])
    VMat[4,:] = np.array([6.7, 0.7, -2.2, -70.7, 12480.0, 81.1, -1.3])
    VMat[5,:] = np.array([-13.7, 11.8, -9.6, -17.0, 81.1, 12630.0, 39.7])
    VMat[6,:] = np.array([-9.9, 4.3, 6.0, -63.3, -1.3, 39.7, 12440.0])
    # convert to a.u.
    VMat = VMat * cm_to_au

    for j in range(NStates):
        VMat[j,j] += np.sum(c[j, :] * R[j*F : (j+1)*F])

    return VMat

def dHel0(R):

    NStates = parameters.NStates
    ndof = parameters.ndof
    F = parameters.F
    ω = parameters.ω
    dVMat0 = np.zeros(( ndof ))

    for j in range(NStates):
        dVMat0[j*F : (j+1)*F] = ω[j, :]**2 * R[j*F : (j+1)*F]

    return dVMat0

def dHel(R):
    
    c = parameters.c
    F = parameters.F
    NStates = parameters.NStates
    ndof    = parameters.ndof

    dVMat = np.zeros(( NStates, NStates, ndof ))

    # diagonal coupling
    for i in range(NStates):
        dVMat[i,i,i*F : (i+1)*F] = c[i,:]

    return dVMat

def initR():

    β  = parameters.β
    ω = parameters.ω
    ndof = parameters.ndof
    F = parameters.F
    NStates = parameters.NStates
    R = np.zeros( ndof )
    P = np.zeros( ndof )

    sigP = np.zeros((NStates, F), dtype=float)
    sigR = np.zeros((NStates, F), dtype=float)

    for i in range(NStates):
        sigP[i, :] = np.sqrt( ω[i, :] / ( 2 * np.tanh( 0.5 * β * ω[i, :] ) ) )
        sigR[i, :] = sigP[i, :] / ω[i, :]     
        for d in range(F):
            R[i*F + d] = np.random.normal() * sigR[i, d]
            P[i*F + d] = np.random.normal() * sigP[i, d]
        
    return R, P