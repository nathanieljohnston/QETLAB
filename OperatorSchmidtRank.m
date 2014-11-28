%%  OPERATORSCHMIDTRANK    Computes the operator Schmidt rank of a bipartite operator
%   This function has one required argument:
%     X: a bipartite operator
%
%   RNK = OperatorSchmidtRank(X) is the operator Schmidt rank of X, which
%   is assumed to live in bipartite space, where both subsystems have
%   dimension equal to sqrt(length(X)).
%
%   This function has one optional argument:
%     DIM (default has two subsystems of equal dimension)
%
%   RNK = OperatorSchmidtRank(X,DIM) is the operator Schmidt rank of X,
%   where the two subsystems it lives on have dimensions specified by the
%   1-by-2 vector DIM.
%
%   URL: http://www.qetlab.com/OperatorSchmidtRank

%   requires: opt_args.m, PermuteSystems.m, SchmidtRank.m, sporth.m, Swap.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: November 12, 2014

function rnk = OperatorSchmidtRank(X,varargin)

dX = size(X);
sdX = round(sqrt(dX));

% set optional argument defaults: dim=sqrt(length(X))
[dim] = opt_args({ [sdX(1) sdX(1);sdX(2) sdX(2)] },varargin{:});

% allow the user to enter a single number for dim
if(length(dim) == 1)
    dim = [dim,dX(1)/dim];
    if abs(dim(2) - round(dim(2))) >= 2*dX(1)*eps
        error('OperatorSchmidtRank:InvalidDim','If DIM is a scalar, X must be square and DIM must evenly divide length(X); please provide the DIM array containing the dimensions of the subsystems.');
    end
    dim(2) = round(dim(2));
end

% allow the user to enter a vector for dim if X is square
if(min(size(dim)) == 1)
    dim = dim(:)'; % force dim to be a row vector
    dim = [dim;dim];
end

% The operator Schmidt rank is just the Schmidt rank of a related vector
% obtained by moving matrix elements around.
rnk = SchmidtRank(Swap(reshape(X,prod(prod(dim)),1),[2,3],[dim(2,:),dim(1,:)]),prod(dim));