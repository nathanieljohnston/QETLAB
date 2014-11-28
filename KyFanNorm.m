%%  KYFANNORM  Computes the Ky Fan k-norm of an operator
%   This function has two required arguments:
%     X: an operator
%     K: a positive integer
%
%   NRM = KyFanNorm(X,K) is the Ky Fan K-norm of X.
%
%   URL: http://www.qetlab.com/KyFanNorm

%   requires: kpNorm.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: December 1, 2012

function nrm = KyFanNorm(X,k)

nrm = kpNorm(X,k,1);