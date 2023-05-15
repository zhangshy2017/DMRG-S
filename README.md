# Extracting Quantum Many-Body Scarred Eigenstates with Matrix Product States

This repository provides the source code to implement the DMRG-S algorithm proposed in [<sup>1</sup>](#refer-anchor-1), which accurately extracts quantum many-body scarred eigenstates. The DMRG-S algorithm can access system sizes far beyond the scope of exact diagonalization and assist analytical studies in discovering exact MPS representations of new scars for generic Hamiltonians. 

The DMRG-S algorithm is implemented based on the ITensor library [<sup>2</sup>](#refer-anchor-2) in Julia programming language. The environment setup requires the installation of the `ITensor.jl` package. The source code consists of three Python scripts:
- `Clifford_circuit_operations.py`
- `Readout_error_correction_functions.py`
- `1D_cluster_state_stabilizer.py`

`Clifford_circuit_operations.py` realizes all the operations of Clifford gates and measuerments on computational bases. The single- and two-qubit depolarization channels are also realized by probabilistic quantum trajectories. We mainly follow the details in [<sup>2</sup>](#refer-anchor-2)

`Readout_error_correction_functions.py` are made up with several functions which exactly or approximately correct the readout errors to accurately compute the entanglement witnesses, such as `Stab_prod_func`, `Stab_prod_exact_correct_func` and `Stab_prod_approx_correct_func`.

`1D_cluster_state_stabilizer.py` is the main processing script, which controls the system size of entangled states, the number of running rounds and correction modes, etc. The whole program can be launched by the terminal command:

`python3 1D_cluster_state_stabilizer.py <Number of entangled qubits> <Correction mode>`

where `Correction mode` has 4 options:
- `Both` for the exact and approximate readout error correction
- `Exact` only for the exact readout error correction
- `Approx` only for the approximate readout error correction
- `None` for no readout error correction

The single-qubit, two-qubit and readout error rates can be changed at the beginning part of `1D_cluster_state_stabilizer.py` (later need to be encapsulated).

The output files are stored in the fold `data` in the form of `.csv` with the naming rules: `1D_cluster_state_n=<Number of entangled qubits>_<0(1) stands for the product of even(odd) stabilizer projectors>mod2_<Correction mode>`.

## References:
<div id="refer-anchor-1"></div>
- [1] S.-Y. Zhang, D. Yuan, T. Iadecola, S. Xu, and D.-L. Deng, Extracting quantum many-body scarred eigenstates with matrix product states (2022), arXiv:2211.05140.

<div id="refer-anchor-2"></div>
- [2] M. Fishman, S. R. White, and E. M. Stoudenmire, The ITensor Software Library for Tensor Network Calculations, SciPost Phys. Codebases , 4 (2022).
