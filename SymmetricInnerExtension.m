%%  SYMMETRICINNEREXTENSION    Determines whether or not an operator has a symmetric inner extension
%   This function has one required argument:
%     X: a positive semidefinite matrix
%
%   EX = SymmetricInnerExtension(X) is either 1 or 0, indicating that X is
%   or is not in the cone defined in reference [1] that approximates the
%   set of separable operators from the inside, based on operators with
%   2-copy positive partial transpose Bosonic symmetric extensions.
%
%   This function has four optional arguments:
%     K (default 2)
%     DIM (default has both subsystems of equal dimension)
%     PPT (default 0)
%     TOL (default eps^(1/4))
%
%   [EX,WIT] = SymmetricInnerExtension(X,K,DIM,PPT,TOL) determines whether
%   or not X is in the cone defined in [1] that approximates the set of
%   separable operators from the insides, based on operators with K-copy
%   Bosonic symmetric extensions, which is positive partial transpose if
%   PPT = 1 and does not have to be positive partial transpose if PPT = 0.
%   If X is in this cone (i.e., EX = 1) then WIT is a symmetric extension
%   of the operator sigma_{AB} from [1], and thus acts as a witness that
%   verifies that EX = 1 is correct. If X is not in this cone (i.e., EX =
%   0) then WIT is an operator with trace(WIT*X) = -1 but trace(WIT*Y) >= 0
%   for all operators in the described cone. Note that WIT may not be an
%   entanglement witness!
%
%   K is the desired number of copies of the second subsystem. DIM is a
%   1-by-2 vector containing the dimensions of the subsystems on which X
%   acts. PPT is a flag (either 0 or 1) indicating whether the desired
%   symmetric extension must have positive partial transpose. TOL is the
%   numerical tolerance used when performing the semidefinite program
%   feasibility check.
%
%   URL: http://www.qetlab.com/SymmetricInnerExtension
%
%   References:
%   [1] M. Navascues, M. Owari, and M. B. Plenio. Complete Criterion for
%       Separability Detection. Physical Review Letters, 103:160404, 2009.

%   requires: cvx (http://cvxr.com/cvx/), IsPPT.m, IsPSD.m, jacobi_poly.m,
%             opt_args.m, PartialTrace.m, PartialTranspose.m,
%             PermutationOperator.m, PermuteSystems.m, sporth.m,
%             SymmetricProjection.m
%             
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: September 22, 2014

function [ex,wit] = SymmetricInnerExtension(X,varargin)

lX = length(X);

% set optional argument defaults: k=2, dim=sqrt(length(X)), ppt=1, tol=eps^(1/4)
[k,dim,ppt,tol] = opt_args({ 2, round(sqrt(lX)), 0, eps^(1/4) },varargin{:});

% allow the user to enter a single number for dim
if(length(dim) == 1)
    dim = [dim,lX/dim];
    if abs(dim(2) - round(dim(2))) >= 2*lX*eps
        error('SymmetricInnerExtension:InvalidDim','If DIM is a scalar, it must evenly divide length(X); please provide the DIM array containing the dimensions of the subsystems.');
    end
    dim(2) = round(dim(2));
end

if(k >= 2)
    if(ppt)
        en = min(1-roots(jacobi_poly(dim(2)-2,mod(k,2),floor(k/2)+1)))*dim(2)/(2*(dim(2)-1));
    end
    
    sdp_dim = [dim(1),dim(2)*ones(1,k)];
    sdp_prod_dim = dim(1)*nchoosek(k+dim(2)-1, dim(2)-1);
    P = kron(speye(dim(1)),SymmetricProjection(dim(2),k,1));

    cvx_begin sdp quiet
        cvx_precision(tol);
        variable rho(sdp_prod_dim,sdp_prod_dim) hermitian
        dual variable W
        subject to
            rho >= 0;
            if(ppt)
                W : (1-en)*PartialTrace(P*rho*P',3:k+1,sdp_dim) + en*kron(PartialTrace(P*rho*P',2:k+1,sdp_dim),eye(dim(2)))/dim(2) == X;
            else
                W : PartialTrace(P*rho*P',3:k+1,sdp_dim)*k/(k+dim(2)) + kron(PartialTrace(P*rho*P',2:k+1,sdp_dim),eye(dim(2)))/(dim(2)+k) == X;                
            end
            if(ppt)
                PartialTranspose(P*rho*P',1:ceil(k/2)+1,sdp_dim) >= 0;
            end
    cvx_end

    % format the output of cvx into a slightly more user-friendly form
    if(strcmpi(cvx_status,'Solved'))
        ex = 1;
        wit = P*rho*P';
    elseif(strcmpi(cvx_status,'Inaccurate/Solved'))
        ex = 1;
        wit = P*rho*P';
        warning('SymmetricInnerExtension:NumericalProblems','Numerical problems encountered by CVX. Consider adjusting the tolerance level TOL and re-running the script.');
    elseif(strcmpi(cvx_status,'Infeasible'))
        ex = 0;
        wit = -W;
    elseif(strcmpi(cvx_status,'Inaccurate/Infeasible'))
        ex = 0;
        wit = -W;
        warning('SymmetricInnerExtension:NumericalProblems','Numerical problems encountered by CVX. Consider adjusting the tolerance level TOL and re-running the script.');
    else
        error('SymmetricInnerExtension:NumericalProblems',strcat('Numerical problems encountered (CVX status: ',cvx_status,'). Please try adjusting the tolerance level TOL.'));
    end
else
    error('SymmetricInnerExtension:InvalidK','K must be an integer >= 2.');
end