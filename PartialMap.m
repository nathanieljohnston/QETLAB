%%  PARTIALMAP    Applies a superoperator to a subsystem of an operator
%   This function has two required arguments:
%     X: a matrix
%     PHI: a superoperator
%
%   PHIX = PartialMap(X,PHI) is the operator (id \otimes PHI)(X). That is,
%   it is the result of applying the superoperator PHI to the second
%   subsystem of X, which is assumed to act on two subsystems of equal
%   dimension (if this is not the case, see the optional arguments below).
%   PHI should be provided either as a Choi matrix, or as a cell with
%   either 1 or 2 columns whose entries are its Kraus operators (see full
%   QETLAB documentation for details).
%
%   This function has two optional arguments:
%     SYS (default 2)
%     DIM (default has two subsystems of equal dimension)
%
%   PHIX = PartialMap(X,PHI,SYS,DIM) is the result of applying the
%   superoperator PHI to the SYS subsystem of X, where the dimensions of
%   the (possibly more than 2) subsystems that X acts on are given in the
%   vector DIM.
%
%   URL: http://www.qetlab.com/PartialMap

%   requires: ApplyMap.m, iden.m, MaxEntangled.m, opt_args.m,
%             PermuteSystems.m, Swap.m
%
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: February 21, 2014

function PhiX = PartialMap(X,Phi,varargin)

dX = size(X);

% Set optional argument defaults: sys=2, dim=sqrt(length(X))
[sys,dim] = opt_args({ 2, round(sqrt(dX)).'*ones(1,2) },varargin{:});

% allow the user to enter a single number for dim
if(length(dim) == 1)
    dim = [dim,dX(1)/dim];
    if abs(dim(2) - round(dim(2))) >= 2*dX(1)*eps
        error('PartialMap:InvalidDim','If DIM is a scalar, X must be square and DIM must evenly divide length(X); please provide the DIM array containing the dimensions of the subsystems.');
    end
    dim(2) = round(dim(2));
end

% allow the user to enter a vector for dim if X is square
if(min(size(dim)) == 1)
    dim = dim(:)'; % force dim to be a row vector
    dim = [dim;dim];
end

% Compute the dimensions of the auxiliary subsystems.
prod_dim_r1 = prod(dim(1,1:sys-1));
prod_dim_c1 = prod(dim(2,1:sys-1));
prod_dim_r2 = prod(dim(1,sys+1:end));
prod_dim_c2 = prod(dim(2,sys+1:end));

if(iscell(Phi)) % the superoperator was provided as a cell of Kraus operators
    % Compute the Kraus operators on the full system.
    sPhi = size(Phi);
    if(sPhi(2) == 1 || (sPhi(1) == 1 && sPhi(2) > 2)) % map is CP
        Phi = Phi(:);
        Phi = cellfun(@(x) kron(kron(speye(prod_dim_r1),x),speye(prod_dim_r2)),Phi,'UniformOutput',false);
    else
        Phi(:,1) = cellfun(@(x) kron(kron(speye(prod_dim_r1),x),speye(prod_dim_r2)),Phi(:,1),'UniformOutput',false);
        Phi(:,2) = cellfun(@(x) kron(kron(speye(prod_dim_c1),x),speye(prod_dim_c2)),Phi(:,2),'UniformOutput',false);
    end

    PhiX = ApplyMap(X,Phi);
else % the superoperator was provided as a Choi matrix
    dPhi = size(Phi);
    dim = [prod_dim_r1,prod_dim_r1,dPhi(1)/dim(1,sys),dim(1,sys),prod_dim_r2,prod_dim_r2;prod_dim_c1,prod_dim_c1,dPhi(2)/dim(2,sys),dim(2,sys),prod_dim_c2,prod_dim_c2];
    psi_r1 = MaxEntangled(prod_dim_r1,1,0);
    psi_c1 = MaxEntangled(prod_dim_c1,1,0);
    psi_r2 = MaxEntangled(prod_dim_r2,1,0);
    psi_c2 = MaxEntangled(prod_dim_c2,1,0);
    
    Phi = PermuteSystems(kron(kron(psi_r1*psi_c1.',Phi),psi_r2*psi_c2.'),[1,3,5,2,4,6],dim);
    PhiX = ApplyMap(X,Phi);
end