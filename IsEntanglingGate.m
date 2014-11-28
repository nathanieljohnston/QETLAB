%%  ISENTANGLINGGATE   Determines if a unitary is an entangling gate
%   This function has one required input argument:
%     U: a unitary matrix
%
%   EG = IsEntanglingGate(U) is either 1 or 0, indicating that U is or is
%   not an entangling quantum gate.
%
%   This function has one optional input argument:
%     DIM (default has two subsystems of equal size)
%
%   [EG,WIT] = IsEntanglingGate(U,DIM) determines that U is or is not an
%   entangling gate, as above. DIM is a vector containing the dimensions
%   of the subsystems that U acts on.
%
%   If EG = 0 then WIT is a struct that contains a decomposition of U that
%   verifies that it is not an entangling gate. More specifically, U can be
%   written as a product of the permutation operator specified by WIT.perm
%   and the elementary tensor with decomposition given by WIT.dec.
%
%   If EG = 1 then WIT is a (sparse) product pure state that is entangled
%   by U. Furthermore, WIT will always have 2 or fewer non-zero entries in
%   this case.
%
%   URL: http://www.qetlab.com/IsEntanglingGate

%   requires: opt_args.m, IsProductOperator.m, IsProductVector.m,
%             perm_inv.m, PermuteSystems.m, SchmidtDecomposition.m
%
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: November 28, 2012

function [eg,wit] = IsEntanglingGate(U,varargin)

dU = size(U);
sdU = round(sqrt(dU));

% set optional argument defaults: dim=sqrt(length(U)), k=0
[dim] = opt_args({ [sdU(1) sdU(1);sdU(2) sdU(2)] },varargin{:});

% allow the user to enter a single number for dim
num_sys = length(dim);
if(num_sys == 1)
    dim = [dim,dU(1)/dim];
    if abs(dim(2) - round(dim(2))) >= 2*dU(1)*eps
        error('IsEntanglingGate:InvalidDim','If DIM is a scalar, U must be square and DIM must evenly divide length(U); please provide the DIM array containing the dimensions of the subsystems.');
    end
    dim(2) = round(dim(2));
    num_sys = 2;
end

% allow the user to enter a vector for dim if U is square
if(min(size(dim)) == 1)
    dim = dim(:)'; % force dim to be a row vector
    dim = [dim;dim];
end

% go through each permutation of the subsystems and check to see if U decomposes into a local unitary with respect to that permutation
p = perms(1:num_sys);
num_sys_fac = factorial(num_sys);
for j = 1:num_sys_fac
    [ipo,wit.dec] = IsProductOperator(PermuteSystems(U,p(j,:),dim,1),[dim(1,p(j,:));dim(2,:)]);
    if(ipo) % the given permutation of U is a product operator, so U is not entangling
        eg = 0;
        wit.perm = perm_inv(p(j,:));
        return
    end
end

% If all of the previous tests failed, the gate is entangling. Find a witness (i.e., a product state that is mapped to an entangled state), if requested.
eg = 1;
if(nargout > 1)
    % I don't like using nested for loops in MATLAB, but I'm not sure of a
    % way around it here. At least we only have 4 nested loops rather than
    % num_sys loops.
    for p = 1:num_sys % subsystem with two non-zero entries
        % Create a cell that contains three identity matrices of useful sizes (that vary with p).
        dp = [prod(dim(2,1:p-1)),dim(2,p),prod(dim(2,p+1:end))];
        Id = {speye(dp(1)) speye(dp(2)) speye(dp(3))};
    
        % Choose which one or two entries should be non-zero.
        ind = [kron((1:dim(2,p)).',ones(1,2));nchoosek(1:dim(2,p),2)];
        indend = nchoosek(dim(2,p)+1,2);
        for q = 1:indend
            for j = 1:dp(1) % subsystems with only one-non-zero entry
                for k = 1:dp(3) % subsystems with only one-non-zero entry
                    wit = kron(kron(Id{1}(:,j),sum(Id{2}(:,ind(q,:)),2)),Id{3}(:,k));
                    if(~IsProductVector(U*wit,dim(1,:)))
                        wit = wit/norm(wit);
                        return
                    end
                end
            end
        end
    end
end