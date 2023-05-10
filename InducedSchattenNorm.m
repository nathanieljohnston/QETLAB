%%  INDUCEDSCHATTENNORM    Computes a lower bound of the induced Schatten p->q norm of a superoperator
%   This function has two required arguments:
%     Phi: a superoperator, specified as either a Choi matrix or a cell
%          array of Kraus operators
%     P: a real number >= 1 or Inf
%
%   NRM = InducedSchattenNorm(PHI,P) is a lower bound of the induced
%   Schatten P-norm of the superoperator PHI. This estimate of the norm is
%   computed via a randomized algorithm, and thus running this function
%   multiple times may produce different lower bounds.
%
%   This function has four optional input arguments:
%     Q: a real number >= 1 or Inf (by default, Q = P)
%     DIM: a 1-by-2 vector containing the input and output dimensions of
%          PHI, in that order (not required if the input and output
%          dimensions are equal, or if PHI is specified as a cell array)
%     TOL: numerical tolerance used to determine when the algorithm stops
%          running (default sqrt(eps))
%     X0: a vector that acts as a starting point for the randomized
%         algorithm (default is randomly-generated)
%
%   This function has one optional output argument:
%     X: the best input matrix that was found (i.e., the matrix that
%        maximizes the Q-Schatten norm of Phi(X) subject to the
%        constraint SchattenNorm(X,P) = 1.
%   
%   [NRM,X] = InducedSchattenNorm(PHI,P,Q,DIM,TOL,V0) is a lower bound of
%   the induced Schatten P->Q-norm of the superoperator PHI. This estimate
%   of the norm is computed via a randomized algorithm, and thus running
%   this function multiple times (with different X0) may produce different
%   lower bounds. Smaller values of TOL give better numerical precision,
%   but increase the running time of the algorithm.
%
%   URL: http://www.qetlab.com/InducedSchattenNorm

%   requires: ApplyMap.m, ChoiMatrix.m, KrausOperators.m, opt_args.m,
%             Realignment.m, SchattenNorm.m, superoperator_dims.m
%             
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: January 11, 2016

function [nrm,X] = InducedSchattenNorm(Phi,p,varargin)

% Start by finding a set of Kraus operators for Phi and computing the
% dimensions that Phi acts on.
if(~iscell(Phi)) % don't alter the Kraus operators -- just wastes time
    Phi = KrausOperators(Phi,varargin{2:end});
end
[da,db] = superoperator_dims(Phi,0,varargin{2:end});

% Set optional argument defaults: q=p, dim, tol=10^-8, X0=-1 (randomly-generated X0)
[q,dim,tol,X0] = opt_args({ p, [da,db], sqrt(eps), -1 },varargin{:});

% Quickly compute in some special cases.
if((p == 2 || strcmpi(p,'fro') == 1) && (q == 2 || strcmpi(q,'fro') == 1))
    [~,S,V] = svd(Realignment(ChoiMatrix(Phi,2),dim).');
    nrm = S(1,1); % the induced Frobenius norm is the norm of the natural representation
    X = reshape(V(:,1),dim(1),dim(1)).'; % the X achieving the maximum is just the original map, viewed in a wonky basis
    return
end

% In all other cases, we iterate to compute the induced Schatten norm.

% If the user specified a starting guess v0, parse it; otherwise randomly
% generate one.
randX0 = 1;
if(max(size(X0)) > 1)
    if(length(X0) ~= dim(1))
        warning('InducedSchattenNorm:DimensionMismatch','The initial matrix X0 must be of size DA-by-DA, where PHI is a map acting on DA-by-DA matrices. Using a randomly-generated initial matrix instead.');
    else
        randX0 = 0;
    end
end
if randX0 % generate a random starting matrix X0, if appropriate
    X = randn(dim(1),dim(1));
    if(~isreal(ChoiMatrix(Phi))) % only add imaginary part to X if Phi is not real (just to make output prettier)
        X = X + 1i*randn(dim(1),dim(1));
    end
else
    X = X0;
end
X = X/SchattenNorm(X,p); % normalize the starting matrix

% Preparation is done; now do the actual iteration.
it_err = 2*tol+1;
Y = ApplyMap(X,Phi);
nrm = SchattenNorm(Y,q);

while it_err > tol
    % First, find the best left matrix Y, keeping the right matrix X fixed.
    [U,S,V] = svd(Y);
    S = diag(S); % only want the diagonal part of S
    
    if(q == Inf)
        [~,ind] = max(S);
        S = zeros(dim(2),1);
        S(ind) = 1;
    else
        S = S/max(S); % pre-process in this way first for numerical reasons
        S = S.^(q-1); % this is the equality condition from the Schatten Holder inequality
        S = S/norm(S,q/(q-1));
    end
    Y = U*diag(S)*V'; % reconstruct the optimal Y from the new SVD that we just computed
    
    % Next, find the best right matrix X, keeping the left marix Y fixed.
    X = ApplyMap(Y,DualMap(Phi,dim));
    [U,S,V] = svd(X);
    S = diag(S); % only want the diagonal part of S
    
    if(p == 1)
        [~,ind] = max(S);
        S = zeros(dim(1),1);
        S(ind) = 1;
    else
        S = S/max(S); % pre-process in this way first for numerical reasons
        S = S.^(1/(p-1)); % this is the equality condition from the Schatten Holder inequality
        S = S/norm(S,p);
    end
    X = U*diag(S)*V'; % reconstruct the optimal X from the new SVD that we just computed
    
    % Check to see if we made any progress; if so, keep iterating.
    Y = ApplyMap(X,Phi);
    new_nrm = SchattenNorm(Y,q);
    it_err = abs(new_nrm - nrm);
    nrm = new_nrm;
end