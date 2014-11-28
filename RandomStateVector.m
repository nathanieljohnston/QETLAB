%%  RANDOMSTATEVECTOR   Generates a random pure state vector
%   This function has one required argument:
%     DIM: the dimension of the Hilbert space that the pure state lives in
%
%   V = RandomStateVector(DIM) generates a DIM-dimensional state vector,
%   uniformly distributed on the (DIM-1)-sphere. Equivalently, these pure
%   states are uniformly distributed according to Haar measure.
%
%   This function has two optional input arguments:
%     RE (default 0)
%     K (default 0)
%
%   V = RandomStateVector(DIM,RE,K) generates a random pure state vector as
%   above. If RE=1 then all coordinates of V will be real. If K=0 then a
%   pure state is generated without considering its Schmidt rank at all. If
%   K>0 then a random bipartite pure state with Schmidt rank <= K is
%   generated (and with probability 1, the Schmidt rank will equal K). If
%   K>0 then DIM is no longer the dimension of the space on which V lives,
%   but rather is the dimension of the *local* systems on which V lives. If
%   these two systems have unequal dimension, you can specify them both by
%   making DIM a 1-by-2 vector containing the two dimensions.
%
%   URL: http://www.qetlab.com/RandomStateVector

%   requires: iden.m MaxEntangled.m, opt_args.m, PermuteSystems.m, Swap.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: November 12, 2014

function v = RandomStateVector(dim,varargin)

% set optional argument defaults: re=0, k=0
[re,k] = opt_args({ 0, 0 },varargin{:});

if(k > 0 && k < min(dim)) % Schmidt rank plays a role
    % allow the user to enter a single number for dim
    if(length(dim) == 1)
        dim = [dim,dim];
    end

    % if you start with a separable state on a larger space and multiply
    % the extra k dimensions by a maximally entangled state, you get a
    % Schmidt rank <= k state
    psi = MaxEntangled(k,1,0);
    a = randn(dim(1)*k,1);
    b = randn(dim(2)*k,1);
    if(~re)
        a = a + 1i*randn(dim(1)*k,1);
        b = b + 1i*randn(dim(2)*k,1);
    end
    v = kron(psi',speye(prod(dim)))*Swap(kron(a,b),[2,3],[k,dim(1),k,dim(2)]);
    v = v/norm(v);
else % Schmidt rank is full, so ignore it
    v = randn(dim,1);
    if(~re)
        v = v + 1i*randn(dim,1);
    end
    v = v/norm(v);
end