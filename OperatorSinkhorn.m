%%  OPERATORSINKHORN    Performs the operator Sinkhorn iteration
%   This function has one required argument:
%     RHO: a density matrix
%
%   SIGMA = OperatorSinkhorn(RHO) is a density matrix that is locally
%   equivalent to RHO, but has both of its partial traces proportional to
%   the identity (RHO must be bipartite; if it is multipartite, see the
%   optional arguments below).
%
%   Such a density matrix SIGMA does not always exist if RHO is low-rank.
%   An error is returned in these cases.
%
%   This function has two optional input arguments:
%     DIM (default has two subsystems of equal dimension)
%     TOL (default sqrt(eps))
%
%   This function has one optional output argument:
%     F: a cell containing local matrices
%
%   [SIGMA,F] = OperatorSinkhorn(RHO,DIM,TOL) returns SIGMA and F such that
%   SIGMA has all of its (single-party) reduced density matrices
%   proportional to the identity, and SIGMA = Tensor(F)*RHO*Tensor(F)'. In
%   other words, F contains invertible local operations that demonstrate
%   that RHO and SIGMA are locally equivalent.
%
%   DIM is a 1-by-2 vector containing the dimensions of the subsystems on
%   which RHO acts (RHO can act on any number of parties). TOL is the
%   numerical tolerance used when determining when the operator Sinkhorn
%   iteration has converged.
%
%   URL: http://www.qetlab.com/OperatorSinkhorn

%   requires: opt_args.m, PartialTrace.m, PermuteSystems.m
%             
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: October 3, 2014

function [sigma,F] = OperatorSinkhorn(rho,varargin)

dX = length(rho);
sdX = round(sqrt(dX));
tr_rho = trace(rho);

% set optional argument defaults: dim=sqrt(length(rho)), tol=sqrt(eps)
[dim,tol] = opt_args({ [sdX, sdX], sqrt(eps) },varargin{:});
num_sys = length(dim);

% allow the user to enter a single number for dim
if(num_sys == 1)
    dim = [dim,dX/dim];
    if abs(dim(2) - round(dim(2))) >= 2*dX*eps
        error('OperatorSinkhorn:InvalidDim','If DIM is a scalar, X must be square and DIM must evenly divide length(X); please provide the DIM array containing the dimensions of the subsystems.');
    end
    dim(2) = round(dim(2));
    num_sys = 2;
end
tr_rho_p = tr_rho^(1/(2*num_sys));

% Prepare the iteration.
for j = num_sys:-1:1
    Prho{j} = eye(dim(j))/dim(j);
    Prho_tmp{j} = Prho{j};
    F{j} = eye(dim(j))*tr_rho_p;
    ldim(j) = prod(dim(1:j-1));
    rdim(j) = prod(dim(j+1:end));
end

% Perform the operator Sinkhorn iteration.
lastwarn(''); % clears any previous warnings
warning('off','MATLAB:singularMatrix'); % we want to catch invertibility warnings, but not display them
warning('off','MATLAB:nearlySingularMatrix'); % we want to catch invertibility warnings, but not display them
it_err = 1;
while it_err > tol
    it_err = 0;
    max_cond = 0;
    
    % Loop over each of the systems and apply a filter on each one.
    try
        for j = 1:num_sys
            % Compute the reduced density matrix on the j-th system.
            Prho_tmp{j} = PartialTrace(rho,setdiff(1:num_sys,j),dim);
            Prho_tmp{j} = (Prho_tmp{j}+Prho_tmp{j}')/2; % for numerical stability
            it_err = it_err + norm(Prho{j}-Prho_tmp{j});
            Prho{j} = Prho_tmp{j};

            % Apply the filter.
            T = sqrtm(inv(Prho{j}))/sqrt(dim(j));
            Tk = kron(speye(ldim(j)),kron(T,speye(rdim(j))));
            rho = Tk*rho*Tk';
            F{j} = T*F{j};

            max_cond = max(max_cond,cond(F{j}));
        end
    catch err
        it_err = 1;
    end
    
    % Make sure that the local transformation performed is invertible --
    % otherwise, the iteration will typically not converge and the results
    % will be meaningless.
    [~,warnid] = lastwarn;
    if(it_err == 1 || max_cond >= 1/tol || strcmpi('MATLAB:nearlySingularMatrix',warnid) || strcmpi('MATLAB:singularMatrix',warnid))
        error('OperatorSinkhorn:LowRank','The operator Sinkhorn iteration does not converge for RHO. This is often the case if RHO is not of full rank.');
    end
end
warning('on','MATLAB:singularMatrix')
warning('on','MATLAB:nearlySingularMatrix')

sigma = (rho+rho')/2; % done for numerical stability reasons
sigma = tr_rho*sigma/trace(sigma); % correct the scaling of the output