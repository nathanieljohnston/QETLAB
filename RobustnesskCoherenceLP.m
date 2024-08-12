%%  ROBUSTNESSKCOHERENCELP    Computes the robustness of k-coherence of a pure quantum state, along with the optimal state sigma via linear programming
%   This function has two required arguments:
%     V: a pure state vector
%     K: a positive integer
%   
%   [ROBK,SIG] = RobustnesskCoherenceLP(V,K) produces the robustness of
%   k-coherence (ROBK) and the closest k-incoherent state (SIG) to a pure
%   state V. SIG is computed via linear programming, which is quicker than
%   the naive method based on semidefinite programming.

%   requires: CVX (http://cvxr.com/cvx/)
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   last updated: May 14, 2018

function [robk,sig] = RobustnesskCoherenceLP(v,k)
    v = sort(abs(v),'descend');
    
    d = length(v);
    rho = v*v';
    [rob,ell] = RobkCohValue(v,k);
    s = sum(v(ell:d));
    p = ell:d;
    pperm = nchoosek(p,k-ell+1);

    cvx_begin quiet
    cvx_precision best;
    variable sig(d,d) symmetric;
    variable tau_coeff(nchoosek(d-ell+1,k-ell+1));
    minimize trace(sig);
    subject to
        diag(sig) >= 0;
        sig - diag(diag(sig)) <= 0;
        sum(sig) == 0;
        
        tau = zeros(d);
        tau_coeff >= 0;
        for j = nchoosek(d-ell+1,k-ell+1):-1:1
            x = [v(1:ell-1);zeros(d-ell+1,1)];
            x(pperm(j,:)) = (s/(k-ell+1))*ones(k-ell+1,1);
            tau = tau + tau_coeff(j)*(x*x');
        end
        rho + sig == tau;
    cvx_end
    
    robk = real(cvx_optval);
    sig = sig/robk;
end