%%  OPERATORSCHMIDTDECOMPOSITION   Computes the operator Schmidt decomposition of a bipartite operator
%   This function has one required argument:
%     X: a bipartite operator
%
%   S = OperatorSchmidtDecomposition(X) is a vector containing the non-zero
%   operator Schmidt coefficients of the bipartite operator X, where the
%   two subsystems are each of size sqrt(length(X)) and X is assumed to be
%   square.
%
%   This function has two optional input arguments:
%     DIM (default has both subsystems of equal size)
%     K (default 0)
%
%   [S,U,V] = OperatorSchmidtDecomposition(X,DIM,K) gives the operator
%   Schmidt coefficients S of the operator X and the corresponding left and
%   right tensor operators in the cells U and V. DIM is a 1x2 vector containing
%   the dimensions of the subsystems that X lives on. If X is non-square,
%   different row and column dimensions can be specified by putting the row
%   dimensions in the first row of DIM and the column dimensions in the 
%   second row of DIM. K is a flag that determines how many terms in the
%   operator Schmidt decomposition should be computed. If K = 0 then all
%   terms with non-zero operator Schmidt coefficients are computed. If
%   K = -1 then all terms (including zero terms) are computed. If K > 0
%   then the K terms with largest operator Schmidt coefficients are
%   computed.
%
%   URL: http://www.qetlab.com/OperatorSchmidtDecomposition

%   requires: opt_args.m, PermuteSystems.m, SchmidtDecomposition.m, Swap.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: November 12, 2014

function [s,U,V] = OperatorSchmidtDecomposition(X,varargin)

dX = size(X);
sdX = round(sqrt(dX));

% set optional argument defaults: dim=sqrt(length(X)), k=0
[dim,k] = opt_args({ [sdX(1) sdX(1);sdX(2) sdX(2)], 0 },varargin{:});

% allow the user to enter a single number for dim
if(length(dim) == 1)
    dim = [dim,dX(1)/dim];
    if abs(dim(2) - round(dim(2))) >= 2*dX(1)*eps
        error('OperatorSchmidtDecomposition:InvalidDim','If DIM is a scalar, X must be square and DIM must evenly divide length(X); please provide the DIM array containing the dimensions of the subsystems.');
    end
    dim(2) = round(dim(2));
end

% allow the user to enter a vector for dim if X is square
if(min(size(dim)) == 1)
    dim = dim(:)'; % force dim to be a row vector
    dim = [dim;dim];
end

% The operator Schmidt decomposition is just the Schmidt decomposition of a
% related vector obtained by moving matrix elements around.
[s,u,v] = SchmidtDecomposition(Swap(reshape(X,prod(prod(dim)),1),[2,3],[dim(2,:),dim(1,:)]),prod(dim),k);

% Now reshape things into the proper output format.
sr = length(s);
U = mat2cell(reshape(u,dim(1,1),dim(2,1)*sr),dim(1,1),dim(2,1)*ones(1,sr));
V = mat2cell(reshape(v,dim(1,2),dim(2,2)*sr),dim(1,2),dim(2,2)*ones(1,sr));