%%  PAULI    Produces a Pauli operator
%   This function has one required argument:
%     IND (an index 0,1,2,3 or 'x', 'y', or 'z')
%
%   P = Pauli(IND) is the 2-by-2 Pauli matrix indicated by the value of
%   IND. IND = 1 gives the Pauli X operator, IND = 2 gives the Pauli Y
%   operator, IND = 3 gives the Pauli Z operator, and IND = 0 gives the
%   identity operator. Alternatively, IND can be set to one of 'x', 'y', or
%   'z' to indicated the Pauli X, Y, or Z operator.
%
%   This function has one optional argument:
%     SP (default 0)
%
%   P = Pauli(IND,SP) is as above, with sparsity of the output determined
%   by the value of SP. If SP = 0 then the output will be full, if SP = 1
%   then the output will be sparse.
%
%   URL: http://www.qetlab.com/Pauli

%   requires: opt_args.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   version: 0.50
%   last updated: December 18, 2013

function p = Pauli(ind,varargin)

% set optional argument defaults: sp=0
[sp] = opt_args({ 0 },varargin{:});

if(ind(1) == 1 || strcmpi(ind,'x'))
    p = [0 1;1 0];
elseif(ind(1) == 2 || strcmpi(ind,'y'))
    p = [0 -1i;1i 0];
elseif(ind(1) == 3 || strcmpi(ind,'z'))
    p = [1 0;0 -1];
else
    p = eye(2);
end

if(sp)
    p = sparse(p);
end