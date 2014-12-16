%%  ISPPT    Determines whether or not a matrix has positive partial transpose
%   This function has one required argument:
%     X: a square matrix
%
%   PPT = IsPPT(X) is either 1 or 0, indicating that X does or does not
%   have positive partial transpose (within numerical error). X is assumed
%   to act on bipartite space.
%
%   This function has three optional input arguments:
%     SYS (default 2)
%     DIM (default sqrt(length(X)))
%     TOL (default sqrt(eps))
%
%   [PPT,WIT] = IsPPT(X,SYS,DIM,TOL) determines whether or not X has
%   positive partial transpose within the tolerance specified by TOL. DIM
%   DIM is a vector containing the dimensions of the subsystems on which X
%   acts, and SYS is a scalar or vector indicating which subsystems the
%   transpose should be applied on. WIT is the eigenvector corresponding to
%   the minimal eigenvalue of the partial transpose of X, and thus can
%   act as a witness that proves X does not have positive partial transpose
%   (i.e., WIT'*PartialTranspose(X,SYS,DIM)*WIT < 0).
%
%   URL: http://www.qetlab.com/IsPPT

%   requires: IsPSD.m, opt_args.m, PartialTranspose.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: December 12, 2014

function [ppt,wit] = IsPPT(X,varargin)

% set optional argument defaults: sys=2, dim=sqrt(length(X)), tol=sqrt(eps)
[sys,dim,tol] = opt_args({ 2, round(sqrt(length(X))), sqrt(eps) },varargin{:});

% Allow this function to be called within CVX optimization problems.
if(isa(X,'cvx'))
    cvx_begin sdp quiet
    subject to
    	PartialTranspose(X,sys,dim) >= 0;
    cvx_end
    ppt = 1-min(cvx_optval,1); % CVX-safe way to map (0,Inf) to (1,0)
    
% If the function is just being called on a non-CVX variable, just check
% the PPT condition normally (which is much faster).
else
    if(nargout > 1)
        [ppt,wit] = IsPSD(PartialTranspose(X,sys,dim),tol);
    else
        ppt = IsPSD(PartialTranspose(X,sys,dim),tol);
    end
end