%%  RelEntCoherence    Computes the relative entropy of coherence of a quantum state
%   This function has one required argument:
%     RHO: a density matrix or pure state vector
%
%   REC = RelEntCoherence(RHO) is the relative entropy of coherence, which
%   equals S(diag(RHO)) - S(RHO), where S(.) is the von Neumann entropy.
%
%   URL: http://www.qetlab.com/RelEntCoherence

%   requires: Entropy.m, pure_to_mixed.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: January 12, 2016

function REC = RelEntCoherence(rho)

rho = pure_to_mixed(rho); % Let the user specify rho as either a density matrix or a pure state vector

REC = Entropy(diag(diag(rho))) - Entropy(rho);