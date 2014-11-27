%%	GHZSTATE		Produces a GHZ state
%	This function has one required input argument:
%       DIM: the dimension of the state.
%
%   GHZ_STATE = GHZState(DIM) returns the Gisin state
%   described in [1].
%
%   For a system of "n" qubits the GHZ state can be written as 
%           |GHZ> = (|0>^{\otimes n} + |1>^{\otimes n})/sqrt(2)
%   
%	References:
%   [1] Going beyond Bell's theorem. D. Greenberger and M. Horne and A. 
%       Zeilinger. E-print: [quant-ph] arXiv:0712.0921. 2007.
%
%	URL: http://www.qetlab.com/GHZState

%	requires: Nothing
%
% 	author: Vincent Russo (vrusso@uwaterloo.ca)
%           Nathaniel Johnston (nathaniel@njohnston.ca)
%	package: QETLAB 
%	last updated: November 27, 2014

function ghz_state = GHZState( dim )

q0 = zeros(dim^2,1);
q0(1) = 1;

q1 = zeros(dim^2,1);
q1(dim^2) = 1;

ghz_state = (q0 + q1)/sqrt(2);

end