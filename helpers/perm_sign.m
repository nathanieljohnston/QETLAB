%%  PERM_SIGN    Computes the sign of a permutation
%   This function has one required argument:
%     PERM: a permutation vector
%
%   SGN = perm_sign(PERM) is the sign (either -1 or 1) of the permutation
%   PERM. In other words, SGN = (-1)^INV, where INV is the number of
%   inversions contained in PERM.
%
%   URL: http://www.qetlab.com/perm_sign

%   requires: nothing
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: November 21, 2012

function sgn = perm_sign(perm)

Id = speye(length(perm));

% At first glance, the following looks like it should be slow (something
% like O(n^3)). However, MATLAB is able to exploit the extreme sparsity of
% Id(:,perm) to make this computation more like O(n).
sgn = det(Id(:,perm));