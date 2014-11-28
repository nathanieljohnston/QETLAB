%%  REALIGNMENT    Computes the realignment of a bipartite operator
%   This function has one required argument:
%     X: a matrix
%
%   RX = Realignment(X) is the realignment of the matrix X, where it is
%   assumed that the number of rows and columns of X are both perfect
%   squares and both subsystems have equal dimension. The realignment is
%   defined by mapping the operator |ij><kl| to |ik><jl| and extending
%   linearly.
%
%   This function has one optional argument:
%     DIM (default has all subsystems of equal dimension)
%
%   RX = Realignment(X,DIM) gives the realignment of the matrix X,
%   where the dimensions of the two subsystems are given by the vector DIM.
%   If X is non-square, different row and column dimensions can be
%   specified by putting the row dimensions in the first row of DIM and the
%   column dimensions in the second row of DIM.
%
%   URL: http://www.qetlab.com/Realignment

%   requires: opt_args.m, PartialTranspose.m, PermuteSystems.m, Swap.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: November 12, 2014

function RX = Realignment(X,varargin)

dX = size(X);
round_dim = round(sqrt(dX));

% set optional argument defaults: dim = sqrt(length(dim))
[dim] = opt_args({ round_dim' },varargin{:});

% allow the user to enter a single number for dim
if(length(dim) == 1)
    dim = [dim,dX(1)/dim];
    if abs(dim(2) - round(dim(2))) >= 2*dX(1)*eps
        error('Realignment:InvalidDim','If DIM is a scalar, X must be square and DIM must evenly divide length(X); please provide the DIM array containing the dimensions of the subsystems.');
    end
    dim(2) = round(dim(2));
end

% allow the user to enter a vector for dim if X is square
if(min(size(dim)) == 1)
    dim = dim(:)'; % force dim to be a row vector
    dim = [dim;dim];
end

% perform the realignment, using the fact that RX = S*PartialTranspose(S*X,1), where S is the swap operator
RX = Swap(PartialTranspose(Swap(X,[1,2],dim,1),1,[dim(1,2),dim(1,1);dim(2,1),dim(2,2)]),[1,2],[dim(2,1),dim(1,1);dim(1,2),dim(2,2)],1);