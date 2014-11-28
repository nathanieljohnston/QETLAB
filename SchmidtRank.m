%%  SCHMIDTRANK    Computes the Schmidt rank of a bipartite vector
%   This function has one required argument:
%     VEC: a bipartite vector to have its Schmidt rank computed
%
%   RNK = SchmidtRank(VEC) is the Schmidt rank of the vector VEC, assumed
%   to live in bipartite space, where both subsystems have dimension equal
%   to sqrt(length(VEC)).
%
%   This function has two optional arguments:
%     DIM (default [sqrt(length(VEC)),sqrt(length(VEC))])
%     TOL (default sqrt(length(VEC))*eps(norm(VEC)))
%
%   RNK = SchmidtRank(VEC,DIM,TOL) is the Schmidt rank of the bipartite
%   vector VEC, where the two subsystems it lives on have dimensions
%   specified by the 1-by-2 vector DIM and the rank is determined as the
%   number of Schmidt coefficients larger than TOL.
%
%   URL: http://www.qetlab.com/SchmidtRank

%   requires: opt_args.m, sporth.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: December 19, 2012

function rnk = SchmidtRank(vec,varargin)

lv = length(vec);
slv = round(sqrt(lv));

% set optional argument defaults: dim=sqrt(length(vec))
[dim,tol] = opt_args({ slv, slv*eps(norm(vec)) },varargin{:});

% allow the user to enter a single number for dim
if(length(dim) == 1)
    dim = [dim,lv/dim];
    if abs(dim(2) - round(dim(2))) >= 2*lv*eps
        error('SchmidtRank:InvalidDim','The value of DIM must evenly divide length(VEC); please provide a DIM array containing the dimensions of the subsystems.');
    end
    dim(2) = round(dim(2));
end

% compute the Schmidt rank
if(issparse(vec))
    [~,rnk] = sporth(reshape(vec,dim(end:-1:1)),tol);
else
    rnk = rank(reshape(vec,dim(end:-1:1)),tol);
end