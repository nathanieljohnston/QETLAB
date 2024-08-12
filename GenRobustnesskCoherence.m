%%  GENROBUSTNESSKCOHERENCE    Computes the generalized robustness of k-coherence of a mixed quantum state
%   This function has two required arguments:
%     RHO: a mixed quantum state
%     K: a positive integer
%   
%   [ROBK,SIG] = GenRobustnesskCoherence(RHO,K) produces the generalized
%   robustness of k-coherence (ROBK) and the closest state (SIG) to a mixed
%   state RHO. This computation is implemented via semidefinite
%   programming.

%   requires: CVX (http://cvxr.com/cvx/), IskCoherent.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   last updated: May 14, 2018

function [robk,sig] = GenRobustnesskCoherence(rho,k)
    d = size(rho,1);

    cvx_begin sdp quiet
    cvx_precision best;
    variable sig(d,d) hermitian;
    minimize trace(sig)
    subject to
        IskCoherent(rho + sig,k) == 1;
        sig >= 0;
    cvx_end
    
    robk = real(cvx_optval);
    sig = sig/robk;
end