%%  DEPHASINGCHANNEL    Produces a dephasing channel
%   This function has one required argument:
%     DIM: the dimensionality on which the channel acts
%
%   DELTA = DephasingChannel(DIM) is the Choi matrix of the completely
%   dephasing channel that acts on DIM-by-DIM matrices.
%
%   This function has one optional argument:
%     P (default 0)
%   
%   DELTA = DephasingChannel(DIM,P) produces the partially dephasing
%   channel (1-P)*D + P*ID, where D is the completely dephasing channel
%   and ID is the identity channel.
%
%   URL: http://www.qetlab.com/DephasingChannel

%   requires: iden.m, MaxEntangled.m, opt_args.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   last updated: January 22, 2015

function delta = DephasingChannel(dim,varargin)

% set optional argument defaults: p=0
[p] = opt_args({ 0 },varargin{:});

% compute the Choi matrix of the depolarizing channel
psi = MaxEntangled(dim,1,0); % gives a sparse non-normalized state
delta = (1-p)*diag(diag(psi*psi')) + p*(psi*psi');