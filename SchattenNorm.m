%%  SCHATTENNORM  Computes the Schatten p-norm of an operator
%   This function has two required arguments:
%     X: an operator
%     P: a real number >= 1, or Inf
%
%   NRM = SchattenNorm(X,P) is the Schatten P-norm of X.
%
%   URL: http://www.qetlab.com/SchattenNorm

%   requires: kpNorm.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: December 1, 2012

function nrm = SchattenNorm(X,p)

nrm = kpNorm(X,min(size(X)),p);