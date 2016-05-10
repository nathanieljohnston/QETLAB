%%  CONCURRENCE    Computes the concurrence of a 2-qubit state
%   This function has one required argument:
%     RHO: a pure state vector or a density matrix
%
%   C = Concurrence(RHO) is the concurrence of the two-qubit state RHO.
%
%   URL: http://www.qetlab.com/Concurrence

%   requires: Pauli.m, pure_to_mixed.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: February 19, 2016

function c = Concurrence(rho)

rho = pure_to_mixed(rho); % Let the user input either a pure state vector or a density matrix.
[m,n] = size(rho);

if m ~= 4 || n ~= 4
    error('Concurrence:InvalidDim','The concurrence is only defined for two-qubit states (i.e., 4-by-4 density matrices).');
end

% Construct the spin-flipped matrix \tilde{\rho}.
Y = Pauli('Y');
rhotilde = kron(Y,Y)*conj(rho)*kron(Y,Y);

% Compute the concurrence.
lam = sort(abs(sqrt(eig(rho*rhotilde))),'descend');
c = max(lam(1)-lam(2)-lam(3)-lam(4),0);