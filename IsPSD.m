%%  ISPSD    Determines whether or not a matrix is positive semidefinite
%   This function has one required argument:
%     X: a square matrix
%
%   PSD = IsPSD(X) is either 1 or 0, indicating that X is or is not
%   positive semidefinite (within reasonable numerical error).
%
%   This function has one optional input argument:
%     TOL (default eps^(3/4))
%
%   [PSD,WIT] = IsPSD(X,TOL) determines whether or not X is positive
%   semidefinite within the tolerance specified by TOL. WIT is the
%   eigenvector corresponding to the minimal eigenvalue of X, and thus can
%   act as a witness that proves X is not positive semidefinite (i.e.,
%   WIT'*X*WIT < 0).
%
%   URL: http://www.qetlab.com/IsPSD

%   requires: opt_args.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: December 24, 2014

function [psd,wit] = IsPSD(X,varargin)

if(size(X,1) ~= size(X,2))
    psd = 0;
    wit = 0;
    return
end

% set optional argument defaults: tol=eps^(3/4)
[tol] = opt_args({ eps^(3/4) },varargin{:});

% Allow this function to be called within CVX optimization problems.
if(isa(X,'cvx'))
    cvx_begin sdp quiet
    subject to
    	X >= 0;
    cvx_end
    psd = 1-min(cvx_optval,1); % CVX-safe way to map (0,Inf) to (1,0)

% If the function is just being called on a non-CVX variable, just check
% positive semidefiniteness normally (which is much faster).
else
    % only check the Hermitian part of X
    X = (X+X')/2;

    % if the user requested the smallest eigenvector, compute it
    if(nargout > 1)
        if(isreal(X))
            eigs_id = 'SA';
        else
            eigs_id = 'SR';
        end
        [wit,eigval] = eigs(X,1,eigs_id);
        psd = (eigval >= -tol);

    % otherwise, use chol to determine positive semidefiniteness, which is both faster and more accurate
    else
        [~,p] = chol(X+(tol+eps)*speye(length(X)));
        psd = (p == 0);
    end
end