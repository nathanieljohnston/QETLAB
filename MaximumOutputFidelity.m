%%  MAXIMUMOUTPUTFIDELITY    Computes the maximum output fidelity of two quantum channels
%   This function has two required input arguments:
%     PHI,PHI: quantum channels, represented as either Choi matrices or
%              cells of Kraus operators
%
%   MOF = MaximumOutputFidelity(PHI,PSI) is the maximum output fidelity
%   between the two quantum channels PHI and PSI. That is, it is the
%   maximum fidelity between states of the form PHI(RHO) and PSI(SIGMA),
%   where RHO and SIGMA are density matrices.
%
%   This function has two optional input arguments:
%     DIM_PHI,DIM_PSI: 1-by-2 vectors containing the input and output
%                      dimensions of PHI and PSI, respectively
%
%   MOF = MaximumOutputFidelity(PHI,PSI,DIM_PHI,DIM_PSI) is as above, where
%   the input and output dimensions of PHI and PSI are specified in the
%   1-by-2 vectors DIM_PHI and DIM_PSI. DIM_PHI and DIM_PSI should be
%   provided if and only if PHI and PSI are have unequal input and output
%   dimensions and are provided as Choi matrices.
%
%   URL: http://www.qetlab.com/MaximumOutputFidelity

%   requires: ComplementaryMap.m, CVX (http://cvxr.com/cvx/),
%             DiamondNorm.m, KrausOperators.m, superoperator_dims.m
%
%   authors: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: November 24, 2014

function mof = MaximumOutputFidelity(Phi,Psi,varargin)

% Compute the dimensions of PHI and PSI.
[da_phi,db_phi,de_phi] = superoperator_dims(Phi,0,varargin{1:min(end,1)});
[da_psi,db_psi,de_psi] = superoperator_dims(Psi,0,varargin{2:end});

if(da_phi ~= da_psi || db_phi ~= db_psi)
    error('MaximumOutputFidelity:InvalidDims','PHI and PSI must have the same input and output dimensions as each other.');
end
min_de = min(de_phi,de_psi);

% Now construct a new map that we will compute the diamond norm of (this
% diamond norm will be the maximum output fidelity that we seek).
Phi = KrausOperators(Phi,[da_phi,db_phi]);
Psi = KrausOperators(Psi,[da_psi,db_psi]);
new_map = ComplementaryMap([Phi(1:min_de),Psi(1:min_de)]);

% Use the fact that the maximum output fidelity is complementary to the
% diamond norm in a natural way.
mof = DiamondNorm(new_map,[da_phi,de_phi]);