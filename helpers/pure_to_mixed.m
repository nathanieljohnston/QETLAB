%%  pure_to_mixed    Converts a state vector or density matrix representation of a state to a density matrix
%   This function has one required argument:
%     PHI: a density matrix or a pure state vector
%
%   RHO = pure_to_mixed(PHI) is a density matrix representation of PHI,
%   regardless of whether PHI is itself already a density matrix, or if it
%   is a pure state vector.
%
%   URL: http://www.qetlab.com/pure_to_mixed

%   requires: nothing
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: January 12, 2016

function rho = pure_to_mixed(phi)

% Compute the size of phi. If it's already a mixed state, leave it alone.
% If it's a vector (pure state), make it into a density matrix.
[m,n] = size(phi);

if(min(m,n) == 1) % it's a pure state vector
    phi = phi(:);
    rho = phi*phi';
elseif(m == n) % it's a density matrix
    rho = phi;
else % it's neither
    error('pure_to_mixed:InvalidDimensions','PHI must be either a vector or a square matrix.');
end