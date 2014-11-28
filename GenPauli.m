%%  GENPAULI    Produces a generalized Pauli operator
%   This function has three required arguments:
%     IND1 (a nonnegative integer from 0 to DIM-1 inclusive)
%     IND2 (a nonnegative integer from 0 to DIM-1 inclusive)
%     DIM (a positive integer indicating the dimension)
%
%   P = GenPauli(IND1,IND2,DIM) is a DIM-by-DIM unitary operator. More
%   specifically, it is the operator X^IND1*Z^IND2, where X and Z are the
%   "shift" and "clock" operators that naturally generalize the Pauli X and
%   Z operators. These matrices span the entire space of DIM-by-DIM
%   matrices as IND1 and IND2 range from 0 to DIM-1, inclusive.
%
%   This function has one optional argument:
%     SP (default 0)
%
%   P = GenPauli(IND1,IND2,DIM,SP) is as above, with sparsity of the output
%   determined by the value of SP. If SP = 0 then the output will be full,
%   if SP = 1 then the output will be sparse.
%
%   URL: http://www.qetlab.com/GenPauli

%   requires: opt_args.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: December 18, 2013

function p = GenPauli(ind1,ind2,dim,varargin)

% set optional argument defaults: sp=0
[sp] = opt_args({ 0 },varargin{:});

w = exp(2i*pi/dim); % primitive root of unity
X = circshift(speye(dim),1); % shift matrix
Z = spdiags((w.^(0:dim-1)).',0,dim,dim); % clock matrix

p = X^ind1*Z^ind2;

if(~sp)
    p = full(p);
end