%%  RobustnessCoherence    Computes the robustness of coherence of a quantum state
%   This function has one required argument:
%     RHO: a pure state vector or a density matrix
%
%   ROC = RobustnessCoherence(RHO) is the robustness of coherence (as
%   defined in [1,2]) of the quantum state (density matrix) RHO.
%
%   References: [1] C. Napoli, T. R. Bromley, M. Cianciaruso, M. Piani, 
%                   N. Johnston, G. Adesso. Robustness of coherence: An
%                   operational and observable measure of quantum
%                   coherence. Preprint, 2016.
%               [2] C. Napoli, T. R. Bromley, M. Cianciaruso, M. Piani, 
%                   N. Johnston, G. Adesso. Robustness of asymmetry and
%                   coherence of quantum states. Preprint, 2016.
%
%   URL: http://www.qetlab.com/RobustnessCoherence

%   requires: cvx (http://cvxr.com/cvx/), L1NormCoherence.m,
%             pure_to_mixed.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: January 12, 2016

function RoC = RobustnessCoherence(rho)

rho = pure_to_mixed(rho); % Let the user enter either a pure state vector or a density matrix.
n = length(rho);

% If the state is pure or single-qubit, we can compute it faster by
% recalling that it is equal to the l1-norm of coherence.
if(n <= 2 || rank(rho) == 1)
    RoC = L1NormCoherence(rho);
    return
end

% Robustness of coherence is computed by semidefinite programming
cvx_begin sdp quiet
    cvx_precision best;
    variable inc_state(n,1);
    variable sig(n,n) hermitian;
    
    minimize trace(sig);
    
    subject to
        sig >= 0;
        rho + sig == diag(inc_state);
cvx_end

RoC = real(cvx_optval);