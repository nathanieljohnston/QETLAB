%%  INDUCEDMATRIXNORM    Computes a lower bound of the induced p->q norm of a matrix
%   This function has two required arguments:
%     X: a matrix
%     P: a real number >= 1 or Inf
%
%   NRM = InducedMatrixNorm(X,P) is a lower bound of the induced P-norm of
%   the matrix X. This estimate of the norm is computed via a randomized
%   algorithm, and thus running this function multiple times may produce
%   different lower bounds.
%
%   This function has three optional input arguments:
%     Q: a real number >= 1 or Inf (by default, Q = P)
%     TOL: numerical tolerance used to determine when the algorithm stops
%          running (default sqrt(eps))
%     V0: a vector that acts as a starting point for the randomized
%         algorithm (default is randomly-generated)
%
%   This function has one optional output argument:
%     V: the best right-multiplication vector that was found (i.e., the
%        vector that maximizes norm(X*V,Q) subject to norm(V,P) = 1).
%   
%   [NRM,V] = InducedMatrixNorm(X,P,Q,TOL,V0) is a lower bound of the
%   induced P->Q-norm of the matrix X. This estimate of the norm is
%   computed via a randomized algorithm, and thus running this function
%   multiple times (with different V0) may produce different lower bounds.
%   Smaller values of TOL give better numerical precision, but increase the
%   running time of the algorithm.
%
%   URL: http://www.qetlab.com/InducedMatrixNorm

%   requires: opt_args.m
%             
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: January 8, 2016

function [nrm,v] = InducedMatrixNorm(X,p,varargin)

[n,m] = size(X); % size of the matrix

% Set optional argument defaults: q=p, tol=10^-8, v0=-1 (randomly-generated v0)
[q,tol,v0] = opt_args({ p, sqrt(eps), -1 },varargin{:});

% Quickly compute in some special cases.
if(p == 1 && q == 1)
    [nrm,ind] = max(sum(abs(X),1)); % norm is max abs column sum
    v = zeros(m,1);
    v(ind) = 1;
    return
elseif((p == 2 || strcmpi(p,'fro') == 1) && (q == 2 || strcmpi(q,'fro') == 1))
    [~,nrm,v] = svds(X,1); % norm is largest singular value
    return
elseif(p == Inf && q == Inf)
    nrm = max(sum(abs(X),2)); % norm is max abs row sum
    v = ones(m,1);
    return
end

% In all other cases, we iterate to compute the induced matrix norm.

% If the user specified a starting guess v0, parse it; otherwise randomly
% generate one.
randv0 = 1;
if(max(size(v0)) > 1)
    v0 = v0(:); % make sure it's a column vector
    if(length(v0) ~= m)
        warning('InducedMatrixNorm:DimensionMismatch','The initial vector v0 must have length equal to the number of columns of X. Using a randomly-generated intial vector instead.');
    else
        randv0 = 0;
    end
end
if randv0 % generate a random starting vector v0, if appropriate
    v = randn(m,1);
    if(~isreal(X)) % only add imaginary part to v if X is not real (just to make output prettier)
        v = v + 1i*randn(m,1);
    end
else
    v = v0;
end
v = v/norm(v,p); % normalize the starting vector

% Preparation is done; now do the actual iteration.
it_err = 2*tol+1;
nrm = norm(X*v,q);

while it_err > tol
    % First, find the best left vector w, keeping the right vector v fixed.
    w = X*v;
    if(q == Inf)
        [~,ind] = max(abs(w));
        w = zeros(m,1);
        w(ind) = 1;
    else
        wabs = abs(w); % split w into its phases and magnitudes
        wph = w./wabs;
        wph(isnan(wph)) = 1; % take care of division by 0 in previous line
        
        wabs = wabs/max(wabs); % pre-process in this way first for numerical reasons
        wabs = wabs.^(q-1); % this is the equality condition from Holder's inequality
        w = wph.*wabs/norm(wabs,q/(q-1));
    end
    
    % Next, find the best right vector v, keeping the left vector w fixed.
    v = w'*X;
    if(p == 1)
        [~,ind] = max(abs(v));
        v = zeros(n,1);
        v(ind) = 1;
    else
        vabs = abs(v); % split v into its phases and magnitudes
        vph = v./vabs;
        vph(isnan(vph)) = 1; % take care of division by 0 in previous line  
        
        vabs = vabs'/max(vabs); % pre-process in this way first for numerical reasons
        vabs = vabs.^(1/(p-1)); % this is the equality condition from Holder's inequality
        v = vph'.*vabs/norm(vabs,p);
    end
    
    % Check to see if we made any progress; if so, keep iterating.
    new_nrm = norm(X*v,q);
    it_err = abs(new_nrm - nrm);
    nrm = new_nrm;
end