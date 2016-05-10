%%  ENTFORMATION    Computes the entanglement of formation of a bipartite quantum state
%   This function has one required argument:
%     RHO: a density matrix or a pure state vector
%
%   EF = EntFormation(RHO) is the entropy of formation of the bipartite
%   quantum state RHO. Note that this function currently only supports RHO
%   being a pure state or a 2-qubit state: it is not known how to compute
%   the entanglement of formation of higher-dimensional mixed states.
%
%   This function has one optional argument:
%     DIM (default has both subsystems of equal dimension)
%
%   EF = EntFormation(RHO,DIM) is the same as above, where RHO acts on
%   local systems of dimension specified by the 1-by-2 vector DIM.
%
%   URL: http://www.qetlab.com/EntFormation

%   requires: Concurrence.m, Entropy.m, opt_args.m, PartialTrace.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: February 19, 2016

function ef = EntFormation(rho,varargin)

[m,n] = size(rho);
round_dim = round(sqrt(max(m,n)));

% set optional argument defaults: dim = sqrt(length(dim))
[dim] = opt_args({ round_dim' },varargin{:});

% allow the user to enter a single number for dim
if(length(dim) == 1)
    dim = [dim,max(m,n)/dim];
    if abs(dim(2) - round(dim(2))) >= 2*max(m,n)*eps
        error('EntFormation:InvalidDim','If DIM is a scalar, DIM must evenly divide length(RHO); please provide the DIM array containing the dimensions of the subsystems.');
    end
    dim(2) = round(dim(2));
end

if(prod(dim) ~= max(m,n))
    error('EntFormation:InvalidDim','Please provide local dimensions in the argument DIM that match the size of RHO.');
end

% If RHO is a rank-1 density matrix, turn it into a vector instead so we
% can compute the EoF easily.
tmp_rho = orth(rho); % if RHO has rank 1, this matrix has one column, which is the vector RHO projects onto
if(m == n && size(tmp_rho,2) == 1)
    rho = tmp_rho;
    n = 1;
end

% Finally, start actually computing the EoF.
if(min(m,n) == 1) % it's a pure state vector
    rho = rho(:);
    ef = Entropy(PartialTrace(rho*rho',2,dim));
    
elseif(m == n) % it's a density matrix
    if(m == 4) % In the two-qubit case, we know how to compute the EoF exactly.
        % Start by computing the concurrence. The EoF is computed based on it.
        C = Concurrence(rho);

        % Now compute the EoF.
        C1 = (1 + sqrt(1 - C^2))/2;
        C2 = (1 - sqrt(1 - C^2))/2;
        ef = -C1*log2(C1) - C2*log2(C2);
        
    else
        % Hopefully we will replace this error by a numerical method soon.
        error('EntFormation:InvalidDim','We currently only know how to compute the entanglement of formation for two-qubit states (i.e., 4-by-4 density matrices) and pure states.');
    end

else % it's neither
    error('EntFormation:InvalidDim','RHO must be either a vector or a square matrix.');
end