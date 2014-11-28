%%  NEGATIVITY    Computes the negativity of a bipartite density matrix
%   This function has one required argument:
%     RHO: a density matrix
%
%   NEG = Negativity(RHO) is the negativity of the density matrix RHO,
%   assuming that the two subsystems on which RHO acts are of equal
%   dimension (if the local dimensions are unequal, specify them in the
%   optional DIM argument). The negativity of RHO is the sum of the
%   absolute value of the negative eigenvalues of the partial transpose of
%   RHO.
%
%   This function has one optional argument:
%     DIM (default has both subsystems of equal dimension)
%
%   NEG = Negativity(RHO,DIM) is the same as above, where RHO acts on local
%   systems of dimension specified by the 1-by-2 vector DIM.
%
%   URL: http://www.qetlab.com/Negativity

%   requires: kpNorm.m, opt_args.m, PartialTranspose.m, PermuteSystems.m,
%             TraceNorm.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: October 22, 2014

function neg = Negativity(rho,varargin)

dX = size(rho);
round_dim = round(sqrt(dX));

% set optional argument defaults: dim = sqrt(length(dim))
[dim] = opt_args({ round_dim' },varargin{:});

% allow the user to enter a single number for dim
if(length(dim) == 1)
    dim = [dim,dX(1)/dim];
    if abs(dim(2) - round(dim(2))) >= 2*dX(1)*eps
        error('Negativity:InvalidDim','If DIM is a scalar, RHO must be square and DIM must evenly divide length(RHO); please provide the DIM array containing the dimensions of the subsystems.');
    end
    dim(2) = round(dim(2));
end

if(prod(dim) ~= dX(1))
    error('Negativity:InvalidDim','Please provide local dimension in the argument DIM that match the size of RHO.');
end

% Compute the negativity.
neg = (TraceNorm(PartialTranspose(rho,2,dim)) - 1)/2;