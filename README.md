# Extracting Quantum Many-Body Scarred Eigenstates with Matrix Product States

This repository provides the source code to implement the DMRG-S algorithm proposed in [<sup>1</sup>](#refer-anchor-1), which accurately extracts quantum many-body scarred eigenstates. The DMRG-S algorithm can access system sizes far beyond the scope of exact diagonalization and assist analytical studies in discovering exact MPS representations of new scars for generic Hamiltonians. 

The DMRG-S algorithm is implemented based on the ITensor library [<sup>2</sup>](#refer-anchor-2) in Julia programming language. The environment setup requires the installation of the `ITensor.jl` package (for example, ITensors v0.3.24). The source code consists of three Julia scripts:
- `projmpo.jl`
- `abstractprojmpo.jl`
- `dmrgs.jl`

`projmpo.jl` and `abstractprojmpo.jl` include some minor changes compared with the original files in the folder "src/mps/" of the ITensors pakage, which are listed below:

`projmpo.jl` include “product_label” in the struct ProjMPO to overload the "product" function in `abstractprojmpo.jl`;

`abstractprojmpo.jl` include some new methods to contract local tensors from MPO $(H-\xi)^2$ and MPS $|\psi\rangle$ in order to obtain $\mathcal{A}^{\[i,i+1\]}$ and $\tilde{\psi}^{[i,i+1]}$ [<sup>1</sup>](#refer-anchor-1). In addition, the original function "product" is overloaded to implement the opreration for the matrix $\mathcal{A}_{t,\text{eff}}^{[i,i+1]}$ mutiplying a vector;

`dmrgs.jl` contain the main function of DMRG-S and SIMPS method with two-site optimization.

Usage:

1. Bakeup the original files `projmpo.jl` and `abstractprojmpo.jl` in the folder "src/mps/" of the ITensors pakage and replace them with the revised version.
2. Place the file `dmrgs.jl` under the folder "src/mps/"  and add a line 'include("mps/dmrgs.jl")' in the file "src/mps/ITensors.jl".
3. Add methods "dmrgs," and "simps,".


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
