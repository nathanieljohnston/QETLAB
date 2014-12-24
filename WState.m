%%  WSTATE    Generates a (generalized) W-state
%   This function has one required input argument:
%     Q: the number of qubits
%
%   W_STATE = WState(Q) returns the Q-qubit W-state described in [1].
%
%   The W-state on Q qubits is defined by
%       |W> = 1/sqrt(Q)*(|100...0> + |010...0> + ... + |000...1>).
%
%   This function has one optional input argument:
%     COEFF (default [1,1,...,1]/sqrt(Q)): a 1-by-Q vector of coefficients
%
%   W_STATE = WState(Q,COEFF) returns the same state as above, except the
%   coefficient of the j-th term in the sum that defines the W-state is
%   COEFF(j).
%
%   References:
%   [1] Three qubits can be entangled in two inequivalent ways. W. Dur,
%       G. Vidal, and J. I. Cirac. E-print: arXiv:quant-ph/0005115, 2000.
%
%   URL: http://www.qetlab.com/WState

%   requires: opt_args.m
% 	author: Vincent Russo (vrusso@uwaterloo.ca)
%           Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB 
%   last updated: December 24, 2014

function w_state = WState(q,varargin)

% set optional argument defaults: coeff = [1,1,...,1]/sqrt(q)
[coeff] = opt_args({ ones(1,q)/sqrt(q), 1 },varargin{:});

% Do some error checking.
if q < 2
    error('WState:InvalidQ','Q must be at least 2.');
elseif length(coeff) ~= q
    error('WState:InvalidCoeff','COEFF must be a vector of length equal to Q.');
end

w_state = sparse(2^q,1);
for j = 1:q
    w_state(2^(j-1)+1) = coeff(q-j+1);
end