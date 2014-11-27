%%	WSTATE		Produces a W-state
%	This function has one required input argument:
%       DIM: the dimension of the state.
%
%   W_STATE = WState(DIM) returns the W-state described in [1].
%
%   In quantum information, the W-states of n qubits are defined as
%       |W> = 1/sqrt(n) (|100...0> + |010...0> + ... + |00...1>)
%   
%	References:
%   [1] Three qubits can be entangled in two inequivalent ways. W. Dur and
%       G. Vidal and J. I. Cirac. E-print: arXiv:quant-ph/0005115. 2000.
%
%	URL: http://www.qetlab.com/WState

%	requires: Nothing
%
% 	author: Vincent Russo (vrusso@uwaterloo.ca)
%           Nathaniel Johnston (nathaniel@njohnston.ca)
%	package: QETLAB 
%	last updated: November 27, 2014

function w_state = WState( dim )

q = zeros(2^dim,1);
for i=1:dim
    q(2^(i-1)+1) = 1;
end

w_state = 1/sqrt(dim) * q;

end

