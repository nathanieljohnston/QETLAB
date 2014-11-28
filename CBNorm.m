%%  CBNORM    Computes the completely bounded norm of a superoperator
%   This function has one required input argument:
%     PHI: a superoperator
%
%   CB = CBNorm(PHI) is the completely bounded (CB) norm of the
%   superoperator PHI. PHI should be specified as either a cell with one or
%   two columns of Kraus operators, or as a Choi matrix (see online QETLAB
%   tutorial for details about specifying superoperators). If PHI is
%   provided as a Choi matrix with unequal input and output dimensions, a
%   second argument specifying the dimensions should also be provided (see
%   below).
%
%   This function has one optional input argument:
%     DIM (default has both subsystems of equal dimension)
%
%   CB = CBNorm(PHI,DIM) is the CB norm of PHI, as above, where DIM is a
%   1-by-2 vector containing the input and output dimensions of PHI, in
%   that order (equivalently, these are the dimensions of the first and
%   second subsystems of the Choi matrix PHI, in that order).
%
%   URL: http://www.qetlab.com/CBNorm

%   requires: ApplyMap.m, ComplementaryMap.m, ChoiMatrix.m, cvx
%             (http://cvxr.com/cvx/), DiamondNorm.m, DualMap.m, iden.m,
%             IsCP.m, IsHermPreserving.m, IsPSD.m, KrausOperators.m,
%             MaxEntangled.m, opt_args.m, PermuteSystems.m, Swap.m
%
%   author: Nathaniel Johnston (nathaniel@njohnston.ca), based on an
%           algorithm by John Watrous
%   package: QETLAB
%   last updated: November 21, 2014

function cb = CBNorm(Phi,varargin)

if(~iscell(Phi) && ~isa(Phi,'cvx'))
    Phi = KrausOperators(Phi,varargin{:});
end
cb = DiamondNorm(DualMap(Phi));