%%  ROBKCOHVALUE    Computes the robustness of k-coherence of a pure quantum state
%   This function has two required arguments:
%     V: a pure state vector
%     K: a positive integer
%   
%   [ROB,L] = RobkCohValue(V,K) produces the robustness of k-coherence
%   (ROB) as well as the flag L that indicates which "branch" of the
%   associated formula ROB comes from (see Theorem 1 of the associated
%   paper for details).

%   requires: nothing
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   last updated: May 14, 2018

function [rob,l] = RobkCohValue(v,k)
    n = length(v);
    for l = k:-1:2
        s = sum(v(l:n));
        if(v(l-1) >= s/(k-l+1))
            rob = sum(v(1:(l-1)).^2) + s^2/(k-l+1) - 1;
            return
        end
    end
    
    l = 1;
    rob = sum(v)^2/k - 1;
end