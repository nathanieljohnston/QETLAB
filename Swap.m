%%  SWAP    Swaps two subsystems within a state or operator
%   This function has one required argument:
%     X: a vector or matrix to have its subsystems swapped
%
%   SX = Swap(X) swaps the two subsystems of the vector or square matrix X,
%   where it is assumed that both the number of rows of X and the number of
%   columns of X are perfect squares and the subsystems have equal
%   dimension.
%
%   This function has three optional arguments:
%     SYS (default [1,2])
%     DIM (default [sqrt(length(X)),sqrt(length(X))])
%     ROW_ONLY (default 0)
%
%   SX = Swap(X,SYS,DIM,ROW_ONLY) swaps the two subsystems of the vector or
%   matrix X, where the dimensions of the (possibly more than 2) subsystems
%   are given by DIM and the indices of the two subsystems to be swapped
%   are specified in the 1-by-2 vector SYS. If X is non-square and not a
%   vector, different row and column dimensions can be specified by putting
%   the row dimensions in the first row of DIM and the column dimensions in
%   the second row of DIM. If ROW_ONLY is set to 1, then only the rows of X
%   are swapped, but not the columns -- this is equivalent to multiplying X
%   on the left by the corresponding swap operator, but not on the right.
%
%   URL: http://www.qetlab.com/Swap

%   requires: opt_args.m, PermuteSystems.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: November 12, 2014

function SX = Swap(X,varargin)

dX = size(X);
round_dim = round(sqrt(dX));

% set optional argument defaults: sys=[1,2], dim=round(sqrt(length(X))), row_only=0
[sys,dim,row_only] = opt_args({ [1,2], [round_dim(1),round_dim(1);round_dim(2),round_dim(2)], 0 },varargin{:});

% allow the user to enter a single number for dim
num_sys = length(dim);
if(num_sys == 1)
    dim = [dim,dX(1)/dim;dim,dX(2)/dim];
    if abs(dim(1,2) - round(dim(1,2))) + abs(dim(2,2) - round(dim(2,2))) >= 2*prod(dX)*eps
        error('Swap:InvalidDim','The value of DIM must evenly divide the number of rows and columns of X; please provide the DIM array containing the dimensions of the subsystems.');
    end
    dim(1,2) = round(dim(1,2));
    dim(2,2) = round(dim(2,2));
    num_sys = 2;
end

% verify that the input sys makes sense
if any(sys < 1) || any(sys > num_sys)
    error('Swap:InvalidSys','The subsystems in SYS must be between 1 and length(DIM) inclusive.');
elseif(length(sys) ~= 2)
    error('Swap:InvalidSys','SYS must be a vector with exactly two elements.');
end

% swap the indicated subsystems
perm = 1:num_sys;
perm(sys) = perm(sys(end:-1:1));
SX = PermuteSystems(X,perm,dim,row_only);