#First released on 07/08/2024

This is a python-based code for "Generalized discrete truncated Wigner approximation" (GDTWA) approach for nonadiabatic mixed-quantum-classical dynamics. For more details, see the reference below [J. Chem. Phys. 155, 024111 (2021)]:

https://pubs.aip.org/aip/jcp/article/155/2/024111/1064995/Generalized-discrete-truncated-Wigner

To run the code on PC, first modify the "RUNNING.py" file and choose a model desired to test, then type "python3 RUNNING.py" to run it. The output file is the system reduced density matrix. To run the code parallel on HPCC, one can either use OPENMP or some simple scripts to realize parallelization. An example can be found at:

https://github.com/Okita0512/SemiClassical-NAMD

In the "Results" file, the example calculations for the two state pyrazine model, three state benzene model, and 7 state FMO model under 77K are provided. The two state pyrazine model, three state benzene model results perfectly reproduced the results reported in the literature [J. Chem. Phys. 155, 024111 (2021)]. For the 7 state FMO model under 77K, the short-time dynamics outperforms the SpinLSC-W-focused approach, while the long time dynamics suffer from negative population problem (same as SpinLSC due to the effective zero-point energy parameter larger than 0), not satisfying detailed balance relation. 
