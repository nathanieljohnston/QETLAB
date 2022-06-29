%%  NORMALIZE_COLS    Scales the columns of a matrix to have norm 1
%   This function has one required argument:
%     X: a matrix
%
%   Y = normalize_cols(X) is the same as X, except each of its columns has
%   been divided by its length. If a column has length 0, it is left
%   unchanged.
%
%   URL: http://www.qetlab.com/normalize_cols

%   requires: nothing
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: August 23, 2021

function Y = normalize_cols(X)

Y = X./repmat(sqrt(max(sum(abs(X).^2,1),eps)),size(X,1),1);