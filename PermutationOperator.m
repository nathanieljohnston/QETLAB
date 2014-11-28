%%  PERMUTATIONOPERATOR    Produces a unitary operator that permutes subsystems
%   This function has two required arguments:
%     DIM: the dimensions of the subsystems to be permuted
%     PERM: a permutation vector
%
%   P = PermutationOperator(DIM,PERM) is a unitary operator that permutes
%   the order of subsystems according to the permutation vector PERM, where
%   the ith subsystem has dimension DIM(i).
%
%   This function has three optional arguments:
%     INV_PERM (default 0)
%     SP (default 0)
%   
%   P = PermutationOperator(DIM,PERM,INV_PERM,SP) is the same as above, but
%   it implements the inverse permutation of PERM if INV_PERM=1. The
%   permutation operator returned is full if SP=0 and sparse if SP=1.
%
%   URL: http://www.qetlab.com/PermutationOperator

%   requires: iden.m, opt_args.m, PermuteSystems.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: November 28, 2012

function P = PermutationOperator(dim,perm,varargin)

% set optional argument defaults: inv_perm=0, sp=0
[inv_perm, sp] = opt_args({ 0, 0 },varargin{:});

% allow the user to enter a single number for dim
if(length(dim)==1)
    dim = dim*ones(1,max(perm));
end

% swap the rows of Id appropriately
P = PermuteSystems(iden(prod(dim),sp),perm,dim,1,inv_perm);