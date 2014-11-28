%%  IDEN  Computes a sparse or full identity matrix
%   This function has two required arguments:
%     DIM: the number of rows (or columns) of the identity matrix
%     SP: a flag (1 or 0)
%
%   ID = iden(DIM,SP) returns the DIM-by-DIM identity matrix. If SP = 0
%   then ID will be full. If SP = 1 then ID will be sparse.
%
%   Only use this function within other functions to easily get the correct
%   identity matrix. If you always want either the full or the sparse
%   identity matrix, just use MATLAB's built-in eye or speye function.
%
%   URL: http://www.qetlab.com/iden

%   requires: nothing
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: November 28, 2012

function id = iden(dim,sp)

if(sp)
    id = speye(dim);
else
    id = eye(dim);
end