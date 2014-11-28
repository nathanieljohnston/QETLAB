%%  KPNORM    Computes the (k,p)-norm of a vector or matrix
%   This function has three required arguments:
%     X: a vector or matrix
%     K: a positive integer
%     P: a real number >= 1, or Inf
%
%   NRM = kpNorm(X,K,P) is the P-norm of the vector of the K largest
%   elements of the vector X (if X is a vector), or the P-norm of the
%   vector of the K largest singular values of X (if X is a matrix).
%
%   URL: http://www.qetlab.com/kpNorm

%   requires: nothing
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: October 22, 2014

function nrm = kpNorm(X,k,p)

sX = size(X);
nX = min(sX);
xX = max(sX);
    
% If X is a CVX variable, try to compute the norm in a way that won't piss
% off MATLAB.
if(isa(X,'cvx') == 1)
    if(((nX > 1 && k >= nX) || (nX == 1 && k >= xX)) && p == 2) % Frobenius norm or Euclidean norm of a vector
        nrm = norm(X,'fro');
    elseif(nX > 1 && (p == Inf || k == 1)) % operator norm
        nrm = norm(X);
    elseif(nX > 1 && k >= nX && p == 1) % trace norm
        nrm = norm_nuc(X);
    elseif(nX > 1 && p == 1) % Ky Fan k-norm
        nrm = lambda_sum_largest([zeros(sX(1)),X;X',zeros(sX(2))],min(k,nX));
    elseif(nX == 1 && k >= xX) % vector p-norm
        nrm = norm(X,p);
    elseif(nX == 1 && p == Inf) % k is irrelevant here: all p=Inf vector norms are the same
        nrm = norm(X,Inf);
    elseif(nX == 1 && p == 1) % sum of k largest absolute values
        if(isreal(X))
            nrm = sum_largest([X(:);-X(:)],min(k,xX));
        else % is there a simpler way to handle the complex case??
            nrm = lambda_sum_largest([zeros(xX),diag(X);diag(X)',zeros(xX)],min(k,xX));
        end
    else
        if(nX > 1)
            error('kpNorm:InvalidCVX','The kpNorm function is not able to represent the given norm in a form that CVX understands. Possible choices are the operator norm, trace norm, Frobenius norm, or Ky Fan norms.');
        else
            error('kpNorm:InvalidCVX','The kpNorm function is not able to represent the given norm in a form that CVX understands. Possible choices are the vector p-norms or the sum of the absolute value of the k largest entries of the vector.');
        end
    end
    
% If X is *not* a CVX variable, compute the norm in a faster way that
% won't break for people without CVX.
else
    % If the requested norm is the Frobenius norm, compute it using MATLAB's
    % built-in Frobenius norm calculation, which is significantly faster than
    % computing singular values.
    if((nX > 1 && k >= nX) && p == 2)
        nrm = norm(X,'fro');

    % That's the one and only special case though; otherwise just compute the
    % norm from the singular values.
    else
        % p = Inf corresponds to the operator norm, which is calculated faster
        % via k = p = 1.
        if(p == Inf)
            k = 1;
            p = 1;
        end

        % Try to determine which method will be faster for computing k singular
        % values: svds (best if X is sparse and k is small) or svd (best if k
        % is large and/or X is full).
        adj = 20 + 1000*(~issparse(X));

        if(nX == 1)
            s = sort(X,'descend');
            s = s(1:min(k,xX));
        elseif(k <= ceil(nX/adj))
            s = svds(X,k);
        else
            s = svd(full(X));
            s = s(1:min(k,nX));
        end
        nrm = norm(s,p);
    end
end