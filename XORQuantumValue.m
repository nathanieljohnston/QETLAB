%%  XORQUANTUMVALUE    Computes the quantum value of a nonlocal binary XOR game
%   This function has two required input arguments:
%     P: a matrix whose (s,t)-entry gives the probability that the referee
%        will give Alice the value s and Bob the value t
%     F: a binary matrix whose (s,t)-entry indicates the winning choice
%        (either 0 or 1) when Alice and Bob receive values s and t from the
%        referee
%
%   WG = XORQuantumValue(P,F) is the quantum value of the XOR game
%   specified by the probability matrix P and the binary matrix of winning
%   values F. That is, it is the optimal probability that Alice and Bob can
%   win the game if they are allowed to share entanglement, but not allowed
%   to communicate during the game itself.
%
%   URL: http://www.qetlab.com/XORQuantumValue

%   requires: cvx (http://cvxr.com/cvx/)
%
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   version: 0.50
%   last updated: October 16, 2014

function wg = XORQuantumValue(p,f)

[s,t] = size(p);

% do some error checking
tol = eps*s^2*t^2;
if(abs(sum(sum(p)) - 1) > tol)
    error('XORQuantumValue:InvalidP','P must be a probability matrix: its entries must sum to 1.');
elseif(min(min(p)) < -tol)
    error('XORQuantumValue:InvalidP','P must be a probability matrix: its entries must be nonnegative.');
elseif(~all([s,t] == size(f)))
    error('XORQuantumValue:InvalidDims','P and F must be matrices of the same size.');
end

% use semidefinite programming to compute the value of the game
P = p.*(1-2*f);
Q = [zeros(s),P;P',zeros(t)];

cvx_begin sdp quiet
    cvx_precision default;
    variable X(s+t,s+t) hermitian
    maximize trace(Q*X)
    subject to
        diag(X) == 1;
        X >= 0;
cvx_end

% The above SDP actually computes the bias of the game. Convert it to the
% value of the game now.
wg = real(cvx_optval)/4 + 1/2;