%%  PURITY    Computes the purity of a quantum state
%   This function has one required argument:
%     RHO: a density matrix
%
%   GAMMA = Purity(RHO) is the purity of the quantum state RHO (i.e., GAMMA
%   is the quantity trace(RHO^2)).
%
%   URL: http://www.qetlab.com/Purity

%   requires: nothing
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: October 15, 2014

function gamma = Purity(rho)

gamma = real(trace(rho^2)); % "real" gets rid of close-to-0 imaginary part