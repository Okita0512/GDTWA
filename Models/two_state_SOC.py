import numpy as np
from scipy import special

class parameters():

   NSteps = 600     #int(2*10**6)
   NTraj = 1
   dtN = 2
   dtE = dtN/20
   NStates = 2
   M = 1000
   initState = 0
   nskip = 5
   ndof = 2

def Hel(R):

    V = np.zeros((2,2), dtype = complex)
    A = 0.02
    alpha = 0.0
    B = 3.0
    W = 5.0

    V0 = A * (1.0 - alpha * np.exp(- B**2 * R[0]**2))
    theta = np.pi * (special.erf(B * R[0]) + 1.0) / 2.0
    phi = - W * R[1]

    V[0,0] = V0 * np.cos( theta )
    V[1,1] = - V[0,0]
    V[0,1] = V0 * np.sin( theta ) * np.exp( - 1.0j * phi )
    V[1,0] = V0 * np.sin( theta ) * np.exp( 1.0j * phi )

    return V

def dHel0(R):
    return 0

def dHel(R):

    dVij = np.zeros((2,2,2), dtype = complex)
    A = 0.02
    alpha = 0.0
    B = 3.0
    W = 5.0

    V0 = A * (1.0 - alpha * np.exp(- B**2 * R[0]**2))
    dV0 = 2 * A * B**2 * alpha * R[0] * np.exp( - B**2 * R[0]**2)
    theta = np.pi * (special.erf(B * R[0]) + 1.0) / 2.0
    dcostheta = - np.sin( theta ) * np.sqrt(np.pi) * B * np.exp( - B**2 * R[0]**2)
    dsintheta = np.cos( theta ) * np.sqrt(np.pi) * B * np.exp( - B**2 * R[0]**2)
    phi = - W * R[1]
    dexp_iphi = - 1.0j * W * np.exp( 1.0j * phi )
    dexp__iphi = 1.0j * W * np.exp( - 1.0j * phi )

    dVij[0,0,0] = dV0 * np.cos( theta ) + V0 * dcostheta
    dVij[0,0,1] = 0.0

    dVij[1,0,0] = ( dV0 * np.sin( theta ) + V0 * dsintheta ) * np.exp( 1.0j * phi )
    dVij[1,0,1] = V0 * np.sin( theta ) * dexp_iphi

    dVij[0,1,0] = ( dV0 * np.sin( theta ) + V0 * dsintheta ) * np.exp( - 1.0j * phi )
    dVij[0,1,1] = V0 * np.sin( theta ) * dexp__iphi

    dVij[1,1,0] = - dVij[0,0,0]
    dVij[1,1,1] = - dVij[0,0,1]

    return dVij

def initR():

    ndof = parameters.ndof

    R0 = np.zeros((ndof))
    P0 = np.zeros((ndof))
    R = np.zeros((ndof))
    P = np.zeros((ndof))
    sigR = np.zeros((ndof))
    sigP = np.zeros((ndof))

    R0[0] = - 3.0
    R0[1] = - 3.0
    P0[0] = 6.0
    P0[1] = 6.0
    sigR[0] = 0.5
    sigR[1] = 0.5
    sigP[0] = 1.0
    sigP[1] = 1.0

    for i in range(ndof):
        R[i] = np.random.normal() * sigR[i] + R0[i]
        P[i] = np.random.normal() * sigP[i] + P0[i]

    return R, P