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

%   requires: IsPSD.m, opt_args.m, PartialTranspose.m, PermuteSystems.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   version: 0.50
%   last updated: November 20, 2012

function [ppt,wit] = IsPPT(X,varargin)

% set optional argument defaults: sys=2, dim=sqrt(length(X)), tol=sqrt(eps)
[sys,dim,tol] = opt_args({ 2, round(sqrt(length(X))), sqrt(eps) },varargin{:});

if(nargout > 1)
    [ppt,wit] = IsPSD(PartialTranspose(X,sys,dim),tol);
else
    ppt = IsPSD(PartialTranspose(X,sys,dim),tol);
end