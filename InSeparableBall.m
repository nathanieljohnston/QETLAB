%%  INSEPARABLEBALL    Checks whether or not an operator is in the ball of separability centered at the maximally-mixed state
%   This function has one required input argument:
%     X: a positive semidefinite matrix, or a vector of the eigenvalues of
%        a positive semidefinite matrix
%
%   ISB = InSeparableBall(X) is either 1 or 0, indicating that X is or is
%   not contained within the ball of separable operators centered at the
%   identity matrix/maximally-mixed state. The size of this ball was
%   derived in [1].
%
%   URL: http://www.qetlab.com/InSeparableBall
%
%   References:
%   [1] L. Gurvits and H. Barnum. Largest separable balls around the
%       maximally mixed bipartite quantum state. Phys. Rev. A, 66:062311,
%       2002. E-print: arXiv:quant-ph/0204159

%   requires: nothing
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: November 14, 2014

function isb = InSeparableBall(X)

sz = size(X);
len = max(sz);

% If X is a vector, turn it into a matrix (we could instead turn every
% matrix into a vector of eigenvalues, but that would make the computation
% take O(n^3) time, instead of the current method, which is O(n^2)).
if(min(sz) == 1) % X is a vector of eigenvalues
    X = sparse(diag(X));
end

% If X has trace 0 or less it can't possibly be in the separable ball.
if(trace(X) < len*eps)
    isb = 0;
    return;
end
X = X/trace(X);

% The following check relies on the fact that we scaled X so that
% trace(X) = 1. The following condition is then exactly the Gurvits-Barnum
% condition.
isb = (norm(X/norm(X,'fro')^2 - speye(len),'fro') <= 1);