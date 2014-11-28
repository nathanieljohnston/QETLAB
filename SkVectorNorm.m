%%  SKVECTORNORM  Computes the s(k)-norm of a vector
%   This function has one required argument:
%     VEC: a bipartite vector to have its s(k)-norm computed
%
%   NRM = SkVectorNorm(VEC) is the maximal inner product of the bipartite
%   vector VEC with a separable pure state. This is equal to the largest
%   Schmidt coefficient of VEC.
%
%   This function has two optional input arguments:
%     K (default 1)
%     DIM (default has both subsystems of equal dimension)
%
%   NRM = SkVectorNorm(VEC,K,DIM) is the maximal inner product of VEC with
%   a state that has Schmidt rank <= K. This quantity is equal to the sum
%   of the squares of the singular values of the partial trace of VEC. It
%   is also equal to the Euclidean norm of the vector of VEC's k largest
%   Schmidt coefficients.
%
%   DIM is a 1x2 vector containing the dimensions of the subsystems that
%   VEC lives on. If DIM is a scalar instead of a vector, then it is
%   assumed that the first subsystem of size DIM and the second subsystem
%   of size length(VEC)/DIM.
%
%   URL: http://www.qetlab.com/SkVectorNorm

%   requires: opt_args.m, SchmidtDecomposition.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: December 2, 2012

function nrm = SkVectorNorm(vec,varargin)

lv = length(vec);

% set optional argument defaults: k=1, dim=sqrt(length(vec))
[k,dim] = opt_args({ 1, round(sqrt(lv)) },varargin{:});

% allow the user to enter a single number for dim
if(length(dim) == 1)
    dim = [dim,lv/dim];
    if abs(dim(2) - round(dim(2))) >= 2*lv*eps
        error('SkVectorNorm:InvalidDim','If DIM is a scalar, it must evenly divide length(VEC); please provide a DIM array containing the dimensions of the subsystems.');
    end
    dim(2) = round(dim(2));
end

% It's faster to just compute the norm of VEC directly if that will give
% the correct answer.
if(k >= min(dim))
    nrm = norm(vec);
else
    nrm = norm(SchmidtDecomposition(vec,dim,k));
end