%%  L1NormCoherence    Computes the l1-norm of coherence of a quantum state
%   This function has one required argument:
%     RHO: a density matrix or pure state vector
%
%   L1C = L1NormCoherence(RHO) is the l1-norm of coherence, which is the
%   sum of the absolute values of the off-diagonal entries of the density
%   matrix RHO (in the standard basis).
%
%   URL: http://www.qetlab.com/L1NormCoherence

%   requires: pure_to_mixed.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: January 12, 2016

function L1C = L1NormCoherence(rho)

rho = pure_to_mixed(rho); % Let the user specify rho as either a density matrix or a pure state vector

L1C = sum(sum(abs(rho))) - trace(rho);