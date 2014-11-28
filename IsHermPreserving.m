%%  ISHERMPRESERVING    Determines whether or not a superoperator is Hermiticity preserving
%   This function has one required argument:
%     PHI: a superoperator
%
%   IHP = IsHermPreserving(PHI) is either 1 or 0, indicating that PHI is or
%   is not Hermiticity preserving (within reasonable numerical error).
%
%   This function has one optional input argument:
%     TOL (default eps^(3/4))
%
%   IHP = IsHermPreserving(PHI,TOL) determines whether or not PHI is
%   Hermiticity preserving within the numerical tolerance specified by TOL.
%
%   URL: http://www.qetlab.com/IsHermPreserving

%   requires: ApplyMap.m, ChoiMatrix.m, iden.m, MaxEntangled.m, opt_args.m,
%             PermuteSystems.m, sporth.m, superoperator_dims.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: January 4, 2013

function hp = IsHermPreserving(Phi,varargin)

if(iscell(Phi) && size(Phi,2) == 1)
    hp = 1;
    return
end

% set optional argument defaults: tol=eps^(3/4)
[tol] = opt_args({ eps^(3/4) },varargin{:});

% Phi is Hermiticity-preseerving iff its Choi matrix is Hermitian.
C = ChoiMatrix(Phi);
sC = size(C);
if(sC(1) ~= sC(2))
    hp = 0;
else
    hp = (max(max(abs(C-C'))) <= tol);
end