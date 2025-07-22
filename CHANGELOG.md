# Change Log
All notable changes to QETLAB will be documented in this file.

## [1.0] - 2025-07-22
### Added
- CliqueNumber: Bounds the clique number (i.e., maximum size of a clique) of a graph.
- Concurrence: Computes the concurrence of a 2-qubit state.
- CopositivePolynomial: Creates a homogenous polynomial whose non-negativity is equivalent to copositivity of a given matrix.
- EntangledSubspace: Constructs a basis of a bipartite r-entangled subspace of any dimension.
- EntFormation: Computes the entanglement of formation of a 2-qubit state or a pure state.
- IsCopositive: Determines whether or not a matrix is copositive.
- MatsumotoFidelity: Computes the Matsumoto fidelity of two density matrices.
- ParallelRepetition: Produces the coefficients of a parallel repetition of a nonlocal game, in full probability notation. Replaces the old NonlocalGameValue.m.
- PauliChannel: Fixed a bug where "probability vectors" summing to less than 1 were being accepted as input to the function.
- PolynomialAsMatrix: Creates a compact fully symmetric matrix representation of a polynomial.
- PolynomialOptimize: Bounds the optimal value of a homogeneous polynomial on the unit sphere.
- PolynomialSOS: Bounds the optimal value of a homogeneous polynomial on the unit sphere via the Sum-Of-Squares hierarchy.
- RandomGraph: Generates the adjacency matrix of a random graph.
- RandomPPTState: Generates a random density matrix with positive partial transpose, and optionally low rank.
- helpers/asum_vector: Creates all vectors with binary entries adding to a given value. Used by AntisymmetricProjection.m.
- helpers/asymind: Creates all vectors with strictly increasing permutations of an input vector.
- helpers/cg2fp: Converts a Bell functional or behaviour in Collins-Gisin notation and converts it to full probability notation. Used in BellInequalityMax.m.
- helpers/chshd: Produces the coefficients of the CHSH-d nonlocal game. Can be used as a test case for BellInequalityMax.
- helpers/dec_to_bin: Converts a decimal number to a binary vector. Replaces de2bi from the Communications toolbox.
- helpers/cg2fc: Converts a Bell functional or behaviour in Collins-Gisin notation notation and converts it to full correlator notation. Used in BellInequalityMax.m.
- helpers/cg2fp: Converts a Bell functional or behaviour in Collins-Gisin notation notation and converts it to full probability notation. Used in BellInequalityMax.m.
- helpers/exp2ind: Looks up a monomial's lexicographical index based on a list of exponents.
- helpers/fc2fp: Converts a Bell functional or behaviour in full correlator notation and converts it to full probability notation. Used in BellInequalityMax.m.
- helpers/fc2cg: Converts a Bell functional or behaviour in full correlator notation and converts it to Collins-Gisin notation. Used in BellInequalityMax.m.
- helpers/ffl: Produces the coefficients of the Fortnow-Feige-Lov�sz nonlocal game. Can be used as a test case for BellInequalityMax.
- helpers/fp2fc: Converts a Bell functional or behaviour in full probability notation and converts it to full correlator notation. Used in BellInequalityMax.m.
- helpers/fp2cg: Converts a Bell functional or behaviour in full probability notation and converts it to Collins-Gisin notation. Used in BellInequalityMax.m.
- helpers/has_band_k_ordering: Determines whether a matrix has bandwidth ≤ k up to symmetric permutation. Used as a new helper check in IskIncoherent.m.
- helpers/glob_ind: Creates a global index from a vector of local indices. Used to be bundled inside of SymmetricProjection.m.
- helpers/pad_array: Pads an array with zeroes. Replaces padarray from the Image Processing toolbox.
- helpers/poly_rand_input: Evaluates a homogeneous polynomial on a random input from the unit sphere.
- helpers/sum_vector: Creates all vectors with non-negative integer entries adding to a given value. Used to be bundled inside of SymmetricProjection.m.
- helpers/symind: Creates all vectors with non-increasing permutations of an input vector.
- helpers/symindfind: Finds the row index of a vector in symind.

### Changed
- AntisymmetricProjection: Dramatically increased speed when using MODE = 0. Changed and standardized the order of the columns when using PARTIAL = 1 and MODE = 0.
- BellInequalityMax: Changed to allow input in full probability, full correlation, or Collins-Gisin notation. See new documentation page for details.
- DiamondNorm: Changed the SDP used in the calculation. This function is now more numerically robust, at the expense of being slightly slower.
- DickeState: Changed the code so that it no longer depends on helpers/unique_perms.
- Entropy: Improved numerical stability so that it no longer frequently returns NaN output.
- GHZState: Now accepts DIM = 1 and/or Q = 1 as input.
- IsBlockPositive: Fixed a numerical tolerance error that would sometimes cause incorrect results to be reported.
- IskIncoherent: Fixed bug with nested CVX optimization and added a bandwidth check to sometimes return early. Also fixed documentation.
- IsSeparable: Fixed a numerical tolerance error that would sometimes cause incorrect results to be reported.
- Negativity: Users can now input either a pure state vector or a density matrix (previously, only density matrices were accepted).
- NonlocalGameValue: Now computes classical value of a game quicker, via algorithm of arXiv:2005.13418
- PartialTrace: Now allows pure state vectors as input, and computes their partial traces (i.e., reduced density matrices) much more quickly.
- PartialTranspose: Fixed bug when partial transposing non-numerical non-square matrices.
- PermuteSystems: Modified so that no computation is performed (so this function is slightly faster) if the permutation requested is the identity.
- SymmetricExtension: If PPT = 1 is specified, now all possible PPT constraints are used instead of just one PPT constraint like in the past.
- SymmetricProjection: Dramatically increased speed when using MODE = 0. Changed and standardized the order of the columns when using PARTIAL = 1 and MODE = 0.
- Tensor: Now works properly (i.e., returns the scalar 1) when M = 0.
- UPB: Now supports the GenTiles1 and GenTiles2 UPBs

### Removed
- NonlocalGameValue: Made redundant by newly-added ParallelRepetition function.
- helpers/unique_perms: Changed to DickeState and SymmetricProjection made this function no longer necessary.

## [0.9] - 2016-01-12
### Added
- BCSGameLB: Computes a lower bound on the quantum value of a binary constraint system (BCS) game.
- BCSGameValue: Computes the maximum value of a binary constraint system (BCS) game. In classical and no-signalling settings, the value computed is exact, but the quantum value is just an upper bound.
- InducedMatrixNorm: Computes a lower bound on the induced p->q norm of a matrix.
- InducedSchattenNorm: Computes a lower bound on the induced Schatten p->q norm of a superoperator.
- L1NormCoherence: Computes the l1-norm of coherence of a quantum state.
- NonlocalGameLB: Computes a lower bound on the quantum value of a two-player non-local game.
- RandomPOVM: Computes a random POVM of a specified size and with a specified number of outcomes.
- RelEntCoherence: Computes the relative entropy of coherence of a quantum state.
- RobustnessCoherence: Computes the robustness of coherence of a quantum state.
- TraceDistanceCoherence: Computes the trace distance of coherence of a quantum state.
- helpers/bcs_to_nonlocal: Converts a description of a binary constraint system (BCS) game into a form that can be presented as a general non-local game.
- helpers/pure_to_mixed: Converts a state vector or density matrix representation of a state to a density matrix.

### Changed
- Entropy: Fixed a bug that would cause NaN output for some low-rank input states.
- kpNorm: Can now be used as the objective function or as a constraint in a CVX optimization problem, regardless of k and p (only certain special values of k and p were supported previously).
- kpNormDual: Can now be used as the objective function or as a constraint in a CVX optimization problem, regardless of k and p (only certain special values of k and p were supported previously).
- NonlocalGameValue: Added the REPT optional input argument, which lets the user specify the number of times that the non-local game will be repeated in parallel.
- SchattenNorm: Can now be used as the objective function or as a constraint in a CVX optimization problem, regardless of p (only p = 1, p = 2, and p = Inf were supported previously).

## [0.8] - 2015-04-13
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