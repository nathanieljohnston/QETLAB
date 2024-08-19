%%  ModTraceDistanceCoherenceMixed    Computes the modified trace distance of coherence of a mixed quantum state
%   This function has one required argument:
%     RHO: a density matrix
%
%   [TDC,D] = ModTraceDistanceCoherenceMixed(X) is the modified trace
%   distance of coherence of the quantum state (density matrix) RHO. The
%   optional output argument D is the optimal (unnormalized) incoherent
%   (diagonal) state in the minimization of ||RHO - D||.

%   requires: YALMIP (https://yalmip.github.io/)
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   last updated: August 21, 2017

function [TDC,D] = ModTraceDistanceCoherenceMixed(rho)

    % Compute the size of RHO.
    n = length(rho);

    % Compute via semidefinite programming.
    options = sdpsettings('verbose',0); % quiet
    
    % Variables.
    d = sdpvar(n,1);
    W1 = sdpvar(n,n,'hermitian','complex');
    W2 = sdpvar(n,n,'hermitian','complex');
    
    % Constraints.
    F = [d >= 0, [W1,rho - diag(d);(rho - diag(d))',W2] >= 0];
    
    % Minimize.
    sol = optimize(F,0.5*(trace(W1)+trace(W2)),options);

    % Interpret results.
    if sol.problem == 0 % Everything went fine. Get solution.
        TDC = TraceNorm(rho - diag(value(d)));
        D = value(d);

    else % Something went wrong (numerical problems?)
        TDC = -1;
    end
end