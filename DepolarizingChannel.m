%%  DEPOLARIZINGCHANNEL    Produces a depolarizing channel
%   This function has one required argument:
%     DIM: the dimensionality on which the channel acts
%
%   DELTA = DepolarizingChannel(DIM) is the Choi matrix of the completely
%   depolarizing channel that acts on DIM-by-DIM matrices.
%
%   This function has one optional argument:
%     P (default 0)
%   
%   DELTA = DepolarizingChannel(DIM,P) produces the partially depolarizing
%   channel (1-P)*D + P*ID, where D is the completely depolarizing channel
%   and ID is the identity channel.
%
%   URL: http://www.qetlab.com/DepolarizingChannel

%   requires: iden.m, MaxEntangled.m, opt_args.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   last updated: March 4, 2014

function delta = DepolarizingChannel(dim,varargin)

% set optional argument defaults: p=0
[p] = opt_args({ 0 },varargin{:});

% compute the Choi matrix of the depolarizing channel
psi = MaxEntangled(dim,1,0); % gives a sparse non-normalized state
delta = (1-p)*speye(dim^2)/dim + p*(psi*psi');