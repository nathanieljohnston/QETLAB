%%  MAXENTANGLED    Produces a maximally entangled bipartite pure state
%   This function has one required argument:
%     DIM: the local dimension
%
%   PSI = MaxEntangled(DIM) is a DIM^2-by-1 column vector representing the
%   standard bipartite maximally entangled pure state \sum_i \ket{ii}
%   (normalized to have norm 1).
%
%   This function has two optional arguments:
%     SP (default 0)
%     NRML (default 1)
%   
%   PSI = MaxEntangled(DIM,SP,NRML) produces a maximally entangled pure
%   state as above that is sparse if SP = 1 and is full if SP = 0. The pure
%   state is normalized to have Euclidean norm 1 if NRML = 1, and it is
%   unnormalized (i.e., each entry in the vector is 0 or 1 and the
%   Euclidean norm of the vector is sqrt(DIM)) if NRML = 0.
%
%   URL: http://www.qetlab.com/MaxEntangled

%   requires: iden.m, opt_args.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: November 28, 2012

function psi = MaxEntangled(dim,varargin)

% set optional argument defaults: sp=0, nrml=1
[sp,nrml] = opt_args({ 0, 1 },varargin{:});

% construct the vector
psi = reshape(iden(dim,sp),dim^2,1);
if(nrml)
    psi = psi/sqrt(dim);
end