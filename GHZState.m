%%	GHZSTATE    Generates a (generalized) GHZ state
%	This function has two required input arguments:
%     DIM: the local dimension
%     Q: the number of parties (qubits/qudits)
%
%   GHZ_STATE = GHZState(DIM,Q) returns Q-partite GHZ state acting on DIM
%   local dimensions, described in [1]. For example, GHZState(2,3) returns
%   the standard 3-qubit GHZ state on qubits. The output of this function
%   is sparse.
%
%   For a system of Q qubits (i.e., DIM = 2), the GHZ state can be
%   written as |GHZ> = (|0>^{\otimes Q} + |1>^{\otimes Q})/sqrt(2).
%
%   This function has one optional input argument:
%     COEFF (default [1,1,...,1]/sqrt(DIM)): a 1-by-DIM vector of coefficients
%
%   GHZ_STATE = GHZState(DIM,Q,COEFF) is as above, but the coefficient of
%   the term |0>^{\otimes Q} is COEFF(1), the coefficient of the term
%   |1>^{\otimes Q} is COEFF(2), and so on.
%
%	References:
%   [1] Going beyond Bell's theorem. D. Greenberger and M. Horne and A. 
%       Zeilinger. E-print: [quant-ph] arXiv:0712.0921. 2007.
%
%	URL: http://www.qetlab.com/GHZState

%	requires: opt_args.m
% 	authors: Vincent Russo (vrusso@uwaterloo.ca)
%            Nathaniel Johnston (nathaniel@njohnston.ca)
%	package: QETLAB 
%	last updated: December 15, 2014

function ghz_state = GHZState(dim,q,varargin)

% set optional argument defaults: coeff = [1,1,..,1]/sqrt(dim)
[coeff] = opt_args({ ones(1,dim)/sqrt(dim) },varargin{:});

% Do some error checking.
if dim < 2
    error('GHZState:InvalidDim','DIM must be at least 2.');
elseif q < 2
    error('GHZState:InvalidQ','Q must be at least 2.');
elseif length(coeff) ~= dim
    error('GHZState:InvalidCoeff','COEFF must be a vector of length equal to DIM.');
end

% Construct the state (and do it in a way that is less memory-intesting
% than naively tensoring things together).
dim_sum = 1;
for j = 1:q-1
    dim_sum = dim_sum + dim^j;
end

ghz_state = sparse(dim^q,1);
for j = 1:dim
    ghz_state((j-1)*dim_sum + 1) = coeff(j);
end