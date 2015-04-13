# Change Log
All notable changes to QETLAB will be documented in this file.

## [0.8] - Changes since 2015-04-13
### Added
- BellInequalityMaxQubits: Approximates the optimal value of a Bell inequality in qubit (i.e., 2-dimensional quantum) settings.
- NonlocalGameValue: Computes the maximum value of a nonlocal game in a classical, quantum, or no-signalling setting.

### Changed
- BellInequalityMax: Bug fix when computing the classical value of a Bell inequality using measurements that have values other than 0, 1, 2, ..., d-1.
- KrausOperators: If the zero map is provided as input, this function now returns a single zero matrix Kraus operator, rather than an empty cell containing no Kraus operators.
- XORGameValue: Bug fix when computing the value of some XOR games with complex entries.

## [0.7] - 2015-01-22
### Added
- BellInequalityMax: Computes the maximum value of a Bell inequality in a classical, quantum, or no-signalling setting.
- BreuerState: Generates a Breuer state, which is a specific family of bound entangled states on even local dimensions.
- DephasingChannel: Produces a (completely or partially) dephasing channel.
- DickeState: Generates a Dicke state.
- GHZState: Generates a (generalized) GHZ state.
- GisinState: Generates a 2-qubit Gisin state.
- HorodeckiState: Generates a Horodecki state, which is a specific family of bound entangled states in (3 \otimes 3)- and (2 \otimes 4)-dimensional spaces.
- LocalDistinguishability: Computes the maximum probability of distinguishing quantum states by means of symmetric-extendible measurements.
- NPAHierarchy: Determines whether or not a set of probabilities satisfy the conditions of the NPA hierarchy, which is a necessary condition for the probabilities to arise from quantum mechanics.
- PauliChannel: Generates a Pauli channel (i.e., a quantum channel with Kraus operators that are multiples of Pauli operators).
- RandomProbabilities: Computes a random probability vector, distributed uniformly on the unit simplex.
- WState: Generates a (generalized) W-state.
- XORGameValue: Computes the classical or quantum value of a nonlocal XOR game (replaces the old functions XORClassicalValue and XORQuantumValue).
- helpers/unique_perms: Computes all distinct permutations of a given vector (the same as unique(perms(V),'rows'), but typically faster and less memory-intensive).
- helpers/update_odometer: Increases the entries of a vector subject to constraints on how large the entries of that vector can be. Useful when you want to have k nested for loops, but k isn't specified beforehand.

### Changed
- AbsPPTConstraints: The DIM input argument is now optional, and the LAM input argument can now either be a vector of eigenvalues or a density matrix (it had to be a vector of eigenvalues before).
- IsAbsPPT: Can now be used directly as a constraint or objective function within other CVX optimization problems.
- IsPPT: Can now be used directly as a constraint or objective function within other CVX optimization problems.
- IsPSD: Can now be used directly as a constraint or objective function within other CVX optimization problems.
- SymmetricExtension: Can now be used directly as a constraint or objective function within other CVX optimization problems.
- SymmetricInnerExtension: Can now be used directly as a constraint or objective function within other CVX optimization problems.
- TensorSum: Is now much faster at tensoring together lots (> 20) of sparse vectors.

### Removed
- XORClassicalValue: Merged into XORGameValue
- XORQuantumValue: Merged into XORGameValue

## [0.6] - 2014-11-27
### Added
- AbsPPTConstraints: Computes the linear matrix inequalities that determine whether or not a mixed state is "absolutely PPT".
- Fidelity: Computes the fidelity of two quantum states, or allows the user to optimize over the fidelity of two quantum states in CVX.
- InSeparableBall: Determines whether or not a given mixed state is within the ball of separability centered at the maximally-mixed state.
- IsAbsPPT: Determines whether or not a mixed state is "absolutely PPT".
- MaximumOutputFidelity: Computes the maximum output fidelity of two quantum channels.
- XORClassicalValue: Computes the classical value of a two-player nonlocal XOR game.
- helpers/superoperator_dims: Computes the input, output, and environment dimensions of a superoperator. Introduced in order to clean up the code in many other functions.

### Changed
- BrauerStates: Reversed the order of the input arguments D and P.
- CBNorm and DiamondNorm: Updated the script so that it can now also be used in the objective function or constraints of other CVX optimization problems (so you can minimize the diamond norm of all channels satisfying a given linear constraint, for example).
- Entropy: Now has a second optional input argument, ALPHA, which allows the user to compute arbitrary Renyi entropies, rather than just the usual von Neumann entropy.
- Pauli: Now allows the user to request a several-qubit Pauli operator, rather than just a 1-qubit Pauli operator. Also, this function now returns output that is sparse by default, rather than full.
- Swap and SwapOperator: Reversed the order of the SYS and DIM optional input arguments.
- Twirl: Can now also perform Pauli twirls.

## [0.5] - 2014-11-06
- Initial release of QETLAB.