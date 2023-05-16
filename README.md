# Extracting Quantum Many-Body Scarred Eigenstates with Matrix Product States

This repository provides the source code to implement the DMRG-S algorithm proposed in [<sup>1</sup>](#refer-anchor-1), which accurately extracts quantum many-body scarred eigenstates (a small fraction of non-thermal excited eigenstates of non-integrable Hamiltomians). The DMRG-S algorithm can access system sizes far beyond the scope of exact diagonalization and assist analytical studies in discovering exact MPS representations of new scars for generic Hamiltonians. 

<div align=center>
<img src="assets/Algorithm.png#pic_center" width='50%'>
</div >
  
The DMRG-S algorithm is implemented based on the ITensor library [<sup>2</sup>](#refer-anchor-2) in Julia programming language. The environment setup requires the installation of the `ITensor.jl` package (with version $\ge$ v0.3.24). The changes to the ITensors library consist of three Julia scripts:
- `projmpo.jl`
- `abstractprojmpo.jl`
- `dmrgs.jl`

`projmpo.jl` and `abstractprojmpo.jl` include some changes to the original files in the folder "src/mps/" of the ITensors pakage, which are listed below:

- `projmpo.jl` includes “product_label” in the struct ProjMPO to overload the `product` function in `abstractprojmpo.jl`;
- `abstractprojmpo.jl` includes some new methods to contract local tensors from MPO $(H-\xi)^2$ and MPS $|\psi\rangle$ in order to obtain $\mathcal{A}^{\[i,i+1\]}$ and $\tilde{\psi}^{[i,i+1]}$ [<sup>1</sup>](#refer-anchor-1). In addition, the original function `product` is overloaded to implement the operation for the matrix $\mathcal{A}_{t,\text{eff}}^{[i,i+1]}$ mutiplying a vector;

`dmrgs.jl` contains the main function of DMRG-S and SIMPS method with two-site DMRG optimization.

## Instructions for setup and usage:

1. Replace the original files `projmpo.jl` and `abstractprojmpo.jl` in the folder `src/mps/` of the ITensors pakage with our revised version.
2. Add the file `dmrgs.jl` into the folder `src/mps/`. Then add a line `include("mps/dmrgs.jl")` in the file `src/ITensors.jl` of the ITensors pakage.
3. Add the methods `dmrgs` and `simps` in the file `src/exports.jl` of the ITensors pakage.
4. Using `julia PXP_dmrgs.jl` in the terminal to run the code.

where `PXP_dmrgs.jl` include several adjustable parameters:
- `initial_energy` for the initial setting of  $\xi_0$
- `psi0` for the initial setting of  $|\psi_0\rangle$
- `N` is the system size
- `maxD` sets the maximum bond dimension
- `minvalue` sets threshold for the variance to update the $\xi_t$

The output files are stored in the fold `data` in the form of `.h5` to store the MPSs during the optimization.

## References:
<div id="refer-anchor-1"></div>
- [1] S.-Y. Zhang, D. Yuan, T. Iadecola, S. Xu, and D.-L. Deng, Extracting quantum many-body scarred eigenstates with matrix product states (2022), arXiv:2211.05140.

<div id="refer-anchor-2"></div>
- [2] M. Fishman, S. R. White, and E. M. Stoudenmire, The ITensor Software Library for Tensor Network Calculations, SciPost Phys. Codebases , 4 (2022).
