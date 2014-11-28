%%  GELLMANN    Produces a Gell-Mann operator
%   This function has one required argument:
%     IND (an integer between 0 and 8, inclusive)
%
%   G = GellMann(IND) is the 3-by-3 Gell-Mann matrix indicated by the value
%   of IND. IND = 0 gives the identity matrix, while values of 1 through 8
%   each indicate one of the other 8 Gell-Mann matrices.
%
%   This function has one optional argument:
%     SP (default 0)
%
%   G = GellMann(IND,SP) is as above, with sparsity of the output
%   determined by the value of SP. If SP = 0 then the output will be full,
%   if SP = 1 then the output will be sparse.
%
%   URL: http://www.qetlab.com/GellMann

%   requires: opt_args.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: December 18, 2013

function g = GellMann(ind,varargin)

% set optional argument defaults: sp=0
[sp] = opt_args({ 0 },varargin{:});

if(ind == 1)
    g = [0 1 0;1 0 0;0 0 0];
elseif(ind == 2)
    g = [0 -1i 0;1i 0 0;0 0 0];
elseif(ind == 3)
    g = [1 0 0;0 -1 0;0 0 0];
elseif(ind == 4)
    g = [0 0 1;0 0 0;1 0 0];
elseif(ind == 5)
    g = [0 0 -1i;0 0 0;1i 0 0];
elseif(ind == 6)
    g = [0 0 0;0 0 1;0 1 0];
elseif(ind == 7)
    g = [0 0 0;0 0 -1i;0 1i 0];
elseif(ind == 8)
    g = [1 0 0;0 1 0;0 0 -2]/sqrt(3);
else
    g = eye(3);
end

if(sp)
    g = sparse(g);
end