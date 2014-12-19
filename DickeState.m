%%  DICKESTATE    Generates a Dicke state
%   This function has one required input argument:
%     N: the number of qubits in the Dicke state
%
%   DICKE_STATE = DickeState(N) returns the N-qubit Dicke state, originally 
%   introduced in [1] for describing light emission from a cloud of atoms. 
%   These symmetric states are in a sense "far from separable". Refer to 
%   [2] for a more recent exposition on Dicke states. The output of this
%   function is sparse.
%
%   This function has three optional input arguments:
%     K (default 1): number of excitations
%     NRML (default 1): a flag, either 1 or 0, specifying normalization
%
%   DICKE_STATE = DickeState(N,K,NRML) is as above, but generates the Dicke
%   state with K-level excitations and is normalized to have norm 1 if
%   NRML = 1 (and has all entries either 0 or 1 if NRML = 0).
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

%	requires: opt_args.m, TensorSum.m, unique_perms.m
% 	authors: Vincent Russo (vrusso@uwaterloo.ca)
%            Nathaniel Johnston (nathaniel@njohnston.ca)
%	package: QETLAB 
%	last updated: December 15, 2014

function dicke_state = DickeState(N,varargin)

% set optional argument defaults: k = 1, nrml = 1
[k,nrml] = opt_args({ 1, 1 },varargin{:});

% Make sure that k is valid.
if k < 0 || k > N
    error('DickeState:InvalidK','K must be an integer between 0 and N, inclusive.');
end

num_terms = nchoosek(N,k);

id = speye(2); % generate the appropriate qubit |0> and |1> states
init_state = [ones(1,k),zeros(1,N-k)];

% Each row of the following matrix corresponds to one of the terms in the
% sum that defined the Dicke state.
perm_matrix = unique_perms(init_state);

% Generate the Dicke state itself from the (0,1)-labelled terms in the
% previous matrix.
for i = N:-1:1
    dicke_terms{i} = id(:,perm_matrix(:,i)+1);
end
dicke_state = TensorSum(dicke_terms{:});

% Normalize the state, if requested.
if(nrml == 1)
    dicke_state = dicke_state / sqrt(num_terms);
end