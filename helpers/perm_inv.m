%%  PERM_INV    Computes the inverse of a permutation
%   This function has one required argument:
%     PERM: a permutation vector
%
%   PI = perm_inv(PERM) is the inverse permutation of PERM.
%
%   URL: http://www.qetlab.com/perm_inv

%   requires: nothing
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: November 28, 2012

function res = perm_inv(perm)

[~,res] = sort(perm);