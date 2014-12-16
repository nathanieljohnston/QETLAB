%%  PAULI    Produces a Pauli operator
%   This function has one required argument:
%     IND (an index 0,1,2,3 or 'I', 'X', 'Y', 'Z', or a vector of such indices)
%
%   P = Pauli(IND) is the 2-by-2 Pauli matrix indicated by the value of
%   IND. IND = 1 gives the Pauli X operator, IND = 2 gives the Pauli Y
%   operator, IND = 3 gives the Pauli Z operator, and IND = 0 gives the
%   identity operator. Alternatively, IND can be set to one of 'I', 'X',
%   'Y', or 'Z' to indicated the Pauli identity, X, Y, or Z operator.
%
%   IND can also be a vector, in which case the output will be a
%   multi-qubit Pauli operator whose action on the k-th qubit is described
%   by the Pauli operator specified in the k-th entry of IND.
%
%   This function has one optional argument:
%     SP (default 1)
%
%   P = Pauli(IND,SP) is as above, with sparsity of the output determined
%   by the value of SP. If SP = 0 then the output will be full, if SP = 1
%   then the output will be sparse.
%
%   URL: http://www.qetlab.com/Pauli

%   requires: opt_args.m, Pauli.m, Tensor.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: November 27, 2014

function p = Pauli(ind,varargin)

% set optional argument defaults: sp=1
[sp] = opt_args({ 1 },varargin{:});

num_qubits = length(ind);

% Did the user just request a one-qubit Pauli operator? OK, make it.
if(num_qubits == 1)
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
    
% Did the user request a multi-qubit Pauli operator? Construct the
% one-qubit ones and then tensor.
else
    for j = num_qubits:-1:1
        p_cell{j} = Pauli(ind(j),sp);
    end
    p = Tensor(p_cell);
end