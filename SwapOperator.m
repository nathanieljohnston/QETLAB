%%  SWAPOPERATOR    Produces a unitary operator that swaps two subsystems
%   This function has one required argument:
%     DIM: the dimensions of the subsystems
%
%   S = SwapOperator(DIM) is the unitary operator that swaps two copies of
%   DIM-dimensional space. If the two subsystems are not of the same
%   dimension, DIM should be a 1-by-2 vector containing the dimension of
%   the subsystems.
%
%   This function has two optional arguments:
%     SYS (default [1,2])
%     SP (default 0)
%   
%   S = SwapOperator(DIM,SYS,SP) is the unitary operator that swaps the two
%   subsystems whose indices are contained in the 1-by-2 vector SYS, where
%   the (possibly more than 2) subsystems have dimension specified by the
%   vector DIM (if just a single number is provided for DIM, then it is
%   assumed that there are max(SYS) subsystems, each of dimension DIM). The
%   swap operator produced is sparse if SP = 1 and is full if SP = 0.
%
%   URL: http://www.qetlab.com/SwapOperator

%   requires: iden.m, opt_args.m, PermuteSystems.m, Swap.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   version: 0.50
%   last updated: November 28, 2012

function S = SwapOperator(dim,varargin)

% set optional argument defaults: sys=[1,2], sp=0
[sys,sp] = opt_args({ [1,2], 0 },varargin{:});

% allow the user to enter a single number for dim
if(length(dim)==1)
    dim = dim*ones(1,max(sys));
end

% swap the rows of Id appropriately
S = Swap(iden(prod(dim),sp),dim,sys,1);