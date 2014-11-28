%%  FILTERNORMALFORM    Computes the filter normal form of an operator
%   This function has one required argument:
%     RHO: a density matrix
%
%   XI = FilterNormalForm(RHO) is a vector of the coefficients in RHO's
%   filter normal form (see Section IV.D of [1]), which are useful for
%   showing that RHO is entangled.
%
%   The filter normal form is not guaranteed to exist if RHO is not full
%   rank. If a filter normal form can not be found, an error is returned.
%
%   This function has two optional input arguments:
%     DIM (default has both subsystems of equal dimension)
%     TOL (default sqrt(eps))
%
%   This function has four optional output arguments:
%     GA,GB: cells of mutually orthonormal matrices
%     FA,FB: invertible matrices
%
%   [XI,GA,GB,FA,FB] = FilterNormalForm(RHO,DIM,TOL) returns XI, GA, GB,
%   FA, FB such that (eye(length(RHO)) + TensorSum(XI,GA,GB))/length(RHO)
%   equals kron(FA,FB)*RHO*kron(FA,FB)'. In other words, FA and FB are
%   matrices implementing the local filter, XI is a vector of coefficients
%   in the filter normal form, and GA and GB are cells of matrices in the
%   tensor-sum decomposition of the filter normal form.
%
%   DIM is a 1-by-2 vector containing the dimensions of the subsystems on
%   which RHO acts. TOL is the numerical tolerance used when constructing
%   the filter normal form.
%
%   URL: http://www.qetlab.com/FilterNormalForm
%
%   References:
%   [1] O. Gittsovich, O. Guehne, P. Hyllus, and J. Eisert. Unifying several
%   separability conditions using the covariance matrix criterion. Phys.
%   Rev. A, 78:052319, 2008. E-print: arXiv:0803.0757 [quant-ph]

%   requires: OperatorSchmidtDecomposition.m, OperatorSinkhorn.m,
%             opt_args.m, PartialTrace.m, PermuteSystems.m,
%             SchmidtDecomposition.m, Swap.m
%             
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: October 3, 2014

function [xi,GA,GB,FA,FB] = FilterNormalForm(rho,varargin)

lrho = length(rho);

% set optional argument defaults: dim=sqrt(length(rho)), tol=sqrt(eps)
[dim,tol] = opt_args({ round(sqrt(lrho)), sqrt(eps) },varargin{:});

% allow the user to enter a single number for dim
if(length(dim) == 1)
    dim = [dim,lrho/dim];
    if abs(dim(2) - round(dim(2))) >= 2*lrho*eps
        error('FilterNormalForm:InvalidDim','If DIM is a scalar, it must evenly divide length(RHO); please provide the DIM array containing the dimensions of the subsystems.');
    end
    dim(2) = round(dim(2));
end

try
    [sigma,F] = OperatorSinkhorn(rho,dim);
catch err
    % Operator Sinkhorn didn't converge.
    if(strcmpi(err.identifier,'OperatorSinkhorn:LowRank'))
        error('FilterNormalForm:NoFNF','The state RHO can not be transformed into a filter normal form. This is often the case if RHO is not of full rank.');
    else
        rethrow(err);
    end
end

% Do some post-processing to make the output more useful and consistent
% with the literature.
pd = prod(dim);
[xi,GA,GB] = OperatorSchmidtDecomposition(sigma - trace(sigma)*eye(pd)/pd,dim);
xi = pd*xi;
FA = F{1};
FB = F{2};