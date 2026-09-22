# An exact-diagonalization program for random Heisenberg-Zeeman spin chains

This C+shell+gnuplot program takes the general Heisenberg-Zeeman
spin-half chain with random exchange couplings (or random Zeeman
magnetic fields), and exactly computes in an arbitrary 
azimuthal spin-sector of the many-body Hilbert space,
certain local spin fluctuations between adjacent eigenstates. The 
measured quantity had been postulated to demonstrate many-body 
localization phenomena in the random-field Heisenberg 
spin chains [1].

The random exchange data showed evidence against localization 
and around the same time, it was theoretically suggested that 
unbroken SU(2) symmetry favors thermalization tendencies [2,3].

The computation works in simple steps:

1. Constructing the large state space associated to a globally chosen
   spin-Z quantum number.

2. Using a Lin-table [4] based indexing algorithm for O(1) look-up of
   the states.

3. Defining the random Heisenberg-Zeeman Hamiltonian in the indexed
   basis and using a dense Lapack diagonalization routine to
   obtain its spectrum and eigenstates.

4. Computing spin-observables across the whole many-body spectrum or
   any chosen part of it to study the tendencies towards spatial
   localization in the eigenstates of the random spin system.


While the majority of these computations is offered by widely
available and better maintained packages today, the 
program's strength lies in its bespoke control at constructing 
the indexable Hilbert subspaces for arbitrary azimuthal spin sectors. 
Extending that dexterity with sparse spectral routines and a
highest-weight construction of the total-spin sectors remains 
the next goal.

## References

[1] A. Pal and D. A. Huse, [Many-body localization phase transition](https://doi.org/10.1103/PhysRevB.82.174411), Phys. Rev. B 82, 174411 (2010).

[2] R. Vasseur, A. C. Potter, and S. A. Parameswaran, [Quantum criticality of hot random spin chains](https://doi.org/10.1103/PhysRevLett.114.217201), Phys. Rev. Lett. 114, 217201 (2015).

[3] A. C. Potter and R. Vasseur, [Symmetry constraints on many-body localization](https://doi.org/10.1103/PhysRevB.94.224206), Phys. Rev. B 94, 224206 (2016).

[4] H. Q. Lin, [Exact diagonalization of quantum-spin models](https://doi.org/10.1103/PhysRevB.42.6561), Phys. Rev. B 42, 6561 (1990).
