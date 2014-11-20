%%	DICKESTATE		Returns the Dicke state
%	This function has two required input argument:
%       N: Number of qubits in Dicke state.
%       K: Number of excitations. 
%
%   DICKE_STATE = DickeState(N,K) returns the Dicke state, originally 
%                 introduced in [1] for describing light emission from a 
%                 cloud of atoms. These symmetric states are in a sense 
%                 "far from separable". Refer to [2] for a more recent 
%                 exposition on Dicke states.
%
%                 Notation for variables is inspired from reference [2].
%
%	References:
%   [1] R.H. Dicke. Coherence in spontaneous radiation processes. Phys. 
%       Rev. 93, 99. 1954.
%       
%   [2] M. Bergmann and O. Guhne. Entanglement criteria for Dicke states.
%       E-print: arXiv:1305.2818 [quant-ph]. 2013.
%   
%   [3] E. Wolfe and S.F. Yelin. Certifying separability in symmetric mixed
%       states, and superradiance. E-print: arXiv:1307.5779 [quant-ph].
%       2013.
%
%	URL: http://www.qetlab.com/DickeState

%	requires: Nothing
%
% 	author: Vincent Russo (vrusso@uwaterloo.ca)
%           Nathaniel Johnston (nathaniel@njohnston.ca)
%	package: QETLAB 
%	version: 
%	last updated:

function [ dicke_state ] = DickeState( N,k )

omega = 1/sqrt(nchoosek(N,k));

q0 = [1;0]; q1 = [0;1];
init_state = [ones(1,k), zeros(1,N-k)];
perm_matrix = unique_perms(init_state);

perms{length(perm_matrix), length(perm_matrix)} = [];
for i=1:length(perm_matrix)
    for j=1:length(perm_matrix(i,:))
        if perm_matrix(i,j) == 0
            perms{i}{j} = q0;
        else
            perms{i}{j} = q1;
        end
    end
end

states{length(perm_matrix)} = [];
for i=1:length(perm_matrix)
   states{i} = TensorSum(perms{i}{:});
end

dicke_state = 0;
for i=1:length(states)
    dicke_state = dicke_state + states{i};
end

dicke_state = omega * dicke_state;

end