import numpy as np

# ================================= 1, choose a method ================================
import Methods.GDTWA as method

# ================================= 2, choose a model =================================
# import Models.morse1 as model
# import Models.morse2 as model
# import Models.morse3 as model
# import Models.spinBoson_Ohmic as model
# import Models.spinBoson_Debye as model
# import Models.tully1 as model
# import Models.tully2 as model
# import Models.tully3 as model
# import Models.spinBoson_1D as model
# import Models.pyrazine as model
# import Models.benzene as model  
import Models.FMO as model

# ================ 3, choose a stype for methods that requires a stype ================
stype = "focused"
par =  model.parameters
par.dHel = model.dHel
par.dHel0 = model.dHel0
par.initR = model.initR
par.Hel   = model.Hel
par.stype = stype

NSteps = model.parameters.NSteps
NTraj = model.parameters.NTraj
NStates = model.parameters.NStates
nskip = model.parameters.nskip
dtN = model.parameters.dtN
runTraj = method.runTraj

rho_ensemble = runTraj(par)

dd = 0 if (NSteps % nskip == 0) else 1

# ===================== Produce the Density Matrix =======================
PiiFile = open("DM.txt","w") 
for t in range(NSteps//nskip + dd):
    PiiFile.write(f"{round(t * nskip * dtN, 3)} \t")        # convert to a.u.
    for i in range(NStates):
        for j in range(NStates):
            PiiFile.write(str(rho_ensemble[i,j,t].real / NTraj) + "\t")
            PiiFile.write(str(rho_ensemble[i,j,t].imag / NTraj) + "\t")
    PiiFile.write("\n")
PiiFile.close()

# ================== Produce the von Neumann Entropy ====================
# SFile = open("S.txt","w") 
# SFile.write(f"{0} \t")
# SFile.write(str(0) + "\t")
# SFile.write("\n")
# for t in range(1, NSteps//nskip + dd):
#     SFile.write(f"{round(t * nskip * dtN,3)} \t")
#     Temp = rho_ensemble[:,:,t] / NTraj
#     E = - np.matmul(Temp,logm(Temp))      # np.log(A) is defined as simple logarithm to each element, which is not correct.
#     S = E.trace()           # use logm and expm in scipy to do it
#     SFile.write(str(S.real) + "\t")
#     SFile.write("\n")
# SFile.close()