%%	BREUERSTATE    Produces a Breuer state of even dimension >= 2
%	This function has two required arguments:
%     DIM: the local dimension (for dim >= 2 and even) 
%     LAMBDA: describes the weight of the singlet component
%
%	BREUER_STATE = BreuerState(DIM,LAMBDA) gives a Breuer bound entangled 
%   state for two qudits of local dimension DIM, with the LAMBDA parameter
%   describing the weight of the singlet component, as described in [1].
%   The Breuer state that is returned is sparse.
%
%	References:
%   [1] H-P. Breuer. Optimal entanglement criterion for mixed quantum 
%       states. E-print: arXiv:quant-ph/0605036, 2006.
%
%	URL: http://www.qetlab.com/BreuerState

%	requires: MaxEntangled.m, SymmetricProjection.m
% 	authors: Vincent Russo (vrusso@uwaterloo.ca)
%            Nathaniel Johnston (nathaniel@njohnston.ca)
%	package: QETLAB 
%	last updated: December 15, 2014

function breuer_state = BreuerState(dim, lambda)

% Make sure that the dimension is even.
if mod(dim,2) == 1 || dim <= 0
    error('BreuerState:InvalidDim','DIM must be an even positive integer.');
end

% Start by generating a specific maximally-entangled state PSI.
V = fliplr(sparse(diag((-1).^mod(1:dim,2))));
psi = kron(speye(dim),V)*MaxEntangled(dim,1);

% Mix the maximally-entangled state PSI with the normalized projector to
% the symmetric subspace.
breuer_state = lambda * (psi*psi') + (1-lambda) * 2*SymmetricProjection(dim)/(dim*(dim+1));