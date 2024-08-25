%%  ModTraceDistanceCoherencePure    Computes the modified trace distance of coherence of a pure quantum state
%   This function has one required argument:
%     X: a pure state vector
%
%   [TDC,P] = ModTraceDistanceCoherencePure(X) is the modified trace
%   distance of coherence of the quantum state (vector) X. The optional
%   output argument P is the optimal value P in the minimization of
%   ||x*x' - p*D||, where D is an incoherent (diagonal) state.

%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: August 16, 2017

function [TDC,p] = ModTraceDistanceCoherencePure(x)

    x = sort(abs(x),'descend')/norm(x);
    
    if(x(1) >= 1/sqrt(2)) % Case (b) of the theorem.
        p = 2*x(1)^2 - 1;
        TDC = 2*x(1)*sqrt(1 - x(1)^2);
    else % Case (a) of the theorem.
        p = 0;
        TDC = 1;
    end
end