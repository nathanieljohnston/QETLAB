%%  TRACENORM  Computes the trace norm of an operator
%   This function has one required argument:
%     X: an operator
%
%   NRM = TraceNorm(X) is the trace norm of X.
%
%   URL: http://www.qetlab.com/TraceNorm

%   requires: kpNorm.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: October 22, 2014

function nrm = TraceNorm(X)

nrm = kpNorm(X,min(size(X)),1);