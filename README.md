# Generation and Verification of Genuine Multipartite Entanglement on a Superconducting Quantum Processor

This program numerically evaluates the predicted value of entanglement witnesses for the one-dimensional linear cluster state, with the method of random Clifford circuit simulations. The error correction method in [<sup>1</sup>](#refer-anchor-1) is adopted to mitigate the readout errors.

The source code consists of three Python scripts:
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
- [1] Correcting detection errors in quantum state engineering through data processing, Chao Shen and Lu-Ming Duan, New Journal of Physics 14,053053 (2012)

<div id="refer-anchor-2"></div>
- [2] Improved simulation of stabilizer circuits, Scott Aaronson and Daniel Gottesman, Phys. Rev. A70, 052328 (2004)
