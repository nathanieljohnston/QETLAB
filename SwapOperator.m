%%  SWAPOPERATOR    Produces a unitary operator that swaps two subsystems
%   This function has one required argument:
%     DIM: the dimensions of the subsystems
%
%   S = SwapOperator(DIM) is the unitary operator that swaps two copies of
%   DIM-dimensional space. If the two subsystems are not of the same
%   dimension, DIM should be a 1-by-2 vector containing the dimension of
%   the subsystems.
%
%   This function has one optional argument:
%     SP (default 0)
%   
%   S = SwapOperator(DIM,SP) is as above, but the swap operator produced is
%   sparse if SP = 1 and is full if SP = 0.
%
%   URL: http://www.qetlab.com/SwapOperator

%   requires: iden.m, opt_args.m, PermuteSystems.m, Swap.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: November 12, 2014

function S = SwapOperator(dim,varargin)

% set optional argument defaults: sp=0
[sp] = opt_args({ 0 },varargin{:});

% allow the user to enter a single number for dim
if(length(dim)==1)
    dim = [dim,dim];
end

% swap the rows of Id appropriately
S = Swap(iden(prod(dim),sp),[1,2],dim,1);