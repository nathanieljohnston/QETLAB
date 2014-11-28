%%  SKOPERATORNORM    Bounds the S(k)-norm of an operator
%   This function has one required input argument:
%     X: a square matrix
%
%   LB = SkOperatorNorm(X) is a lower bound of the S(1)-norm of the
%   operator X, which is assumed to act on two systems of the same
%   dimension. Note that X should be positive semidefinite in almost all
%   cases in which this function is used -- the bounds provided are
%   generally not very good otherwise.
%
%   This function has five optional input arguments:
%     K (default 1)
%     DIM (default has both subsystems of equal dimension)
%     STR (default 2)
%     TARGET (default -1)
%     TOL (default eps^(3/8))
%
%   [LB,LWIT,UB,UWIT] = SkOperatorNorm(X,K,DIM,STR,TARGET,TOL) provides a
%   lower bound (LB) and an upper bound (UB) of the S(K)-norm of the
%   operator X, which acts on two subsystems, whose dimensions are given in
%   the 1-by-2 vector DIM.
%
%   K is the "index" of the norm -- that is, it is the Schmidt rank of the
%   vectors that are multiplying X on the left and right in the definition
%   of the norm.
%
%   STR is the amount of computation that you want to devote to computing
%   the bounds. Higher values of STR correspond to more computation time
%   (and thus better bounds).
%
%   TARGET is a target value that you wish to prove that the norm is above
%   or below. Once the script has proved that LB >= TARGET or UB <= TARGET,
%   the script immediately aborts. This is a time-saving feature that can
%   be avoided by setting TARGET = -1.
%
%   TOL is the numerical tolerance used throughout the script.
%
%   LWIT and UWIT are witnesses that verify that the bounds LB and UB are
%   correct. More specifically, LWIT is a vector with Schmidt rank <= K
%   such that LWIT'*X*LWIT = LB, and UWIT is the optimal positive
%   semidefinite operator Y in the dual semidefinite program (5.2.3) in [3]
%   (see online documentation for more details).
%
%   Most of the lower bounds and basic facts about these norms were derived
%   in [1]. The semidefinite program method for computing upper bounds was
%   derived in [2]. See [3] for a summary of results and more information.
%
%   URL: http://www.qetlab.com/SkOperatorNorm
%
%   References:
%   [1] N. Johnston and D. W. Kribs. A Family of Norms With Applications in
%       Quantum Information Theory. Journal of Mathematical Physics,
%       51:082202, 2010.
%   [2] N. Johnston and D. W. Kribs. A Family of Norms With Applications in
%       Quantum Information Theory II. Quantum Information & Computation,
%       11(1 & 2):104-123, 2011.
%   [3] N. Johnston. Norms and Cones in the Theory of Quantum Entanglement.
%       PhD thesis, University of Guelph, 2012.

%   requires: cvx (http://cvxr.com/cvx/), iden.m, IsPSD.m, kpNorm.m,
%             MaxEntangled.m, normalize_cols.m, opt_args.m, PartialMap.m,
%             PartialTrace.m, PartialTranspose.m, PermuteSystems.m,
%             Realignment.m, SchmidtDecomposition.m, SchmidtRank.m,
%             sk_iterate.m, SkVectorNorm.m, sporth.m, Swap.m
%
%   author: Nathaniel Johnston (nathaniel@njohnston.ca), based on joint
%           work with David W. Kribs
%   package: QETLAB
%   last updated: September 22, 2014

function [lb,lwit,ub,uwit] = SkOperatorNorm(X,varargin)

    X = full(X);
    dX = length(X);
    sdX = round(sqrt(dX));
    lwit = 0;
    uwit = 0;

    % Set optional argument defaults: k=1, dim=sqrt(length(X)), str=2,
    % target=-1, tol=eps^(3/8)
    [k,dim,str,target,tol] = opt_args({ 1, [sdX sdX], 2, -1, eps^(3/8) },varargin{:});
    if(str == -1)
        str = 1/eps; % keep going forever!
    end
    
    % allow the user to enter a single number for dim
    if(length(dim) == 1)
        dim = [dim,dX/dim];
        if abs(dim(2) - round(dim(2))) >= 2*dX*eps
            error('SkOperatorNorm:InvalidDim','If DIM is a scalar, X must be square and DIM must evenly divide length(X); please provide the DIM array containing the dimensions of the subsystems.');
        end
        dim(2) = round(dim(2));
    end

    % some useful, repeatedly-used, values
    prod_dim = prod(dim);
    op_norm = norm(X);
    x_rank = rank(X);
    X = X/op_norm; % rescale X to have unit norm

    % The S(k)-norm is just the operator norm if k is large enough.
    if k >= min(dim)
        lb = op_norm;
        ub = op_norm;
        return
    end

    % If X is rank 1 then the S(k)-norm is easy to compute via Proposition
    % 10 of [1].
    if x_rank == 1
        [U,~,V] = svds(X,1);
        lb = op_norm*SkVectorNorm(U,k,dim)*SkVectorNorm(V,k,dim);
        ub = lb;
        return
    end

    isprojection = 0;
    ishermitian = 0;
    ispositive = 0;
    isnormal = (max(max(abs(X'*X-X*X'))) <= 2*prod_dim*eps);
    if isnormal
        ishermitian = (max(max(abs(X-X'))) <= 100*eps);
        if ishermitian
            X = (X+X')/2; % avoids some numerical problems later
            ispositive = IsPSD(X);
            if ispositive
                isprojection = (max(max(abs(X-X*X))) <= eps*prod_dim^2);
            end
        end
    end
    is_trans_exact = (min(dim) == 2 && max(dim) <= 3);

    % Compute some more simple bounds. We will multiply these by op_norm before
    % we leave this function.
    lb = (k/min(dim));    % comes from Theorem 17 in [1]
    ub = 1; % our most basic upper bound

    % break out of the function if the target value has already been met
    if(MetTarget(lb,ub,op_norm,tol,target))
        lb = op_norm*lb;
        ub = op_norm*ub;
        return
    end
    
    if ~(ispositive && is_trans_exact && k == 1 && str >= 1) % if the exact answer won't be found by SDP, compute bounds via other methods first
        if(ishermitian)
            [eig_vec,eig_val] = eig(X);
            [eig_val,ord] = sort(real(diag(eig_val)),'ascend');
            eig_vec = eig_vec(:,ord);

            % use the lower bound of Proposition 4.14 of [1]
            for r = k:min(dim)
                t_ind = prod(dim) - prod(dim - r);
                lb = max([(k/r)*eig_val(t_ind),lb]);
            end

            % use the lower bound of Theorem 4.2.15 of [3]
            if(k==1)
                lb = max([(trace(X) + sqrt((prod_dim*trace(X^2)-trace(X)^2)/(prod_dim-1)))/prod_dim,lb]);
            end
        end

        if(ispositive && nargout > 2)
            % Use the upper bound of Proposition 15 of [1].
            t_ub = 0;
            for i = 1:prod_dim
                t_ub = t_ub + abs(eig_val(i))*SkVectorNorm(eig_vec(:,i),k,dim)^2;
            end
            ub = min(t_ub,ub);

            % Use the upper bound of Proposition 4.2.11 of [3].
            ub = min(kpNorm(Realignment(X,dim),k^2,2),ub);
        end

        % Use the lower bound of Theorem 4.2.17 of [3].
        if(isprojection)
            lb = max(min([1,k/ceil((dim(1)+dim(2) - sqrt((dim(1)-dim(2))^2 + 4*x_rank - 4))/2)]),lb);
            lb = max((min(dim)-k)*(x_rank+sqrt((prod_dim*x_rank-x_rank^2)/(prod_dim-1)))/(prod_dim*(min(dim)-1)) + (k-1)/(min(dim)-1),lb);
        end

        % break out of the function if the target value has already been met
        if(MetTarget(lb,ub,op_norm,tol,target))
            lb = op_norm*lb;
            ub = op_norm*ub;
            return
        end

        % Use a randomized iterative method to try to improve the lower
        % bound.
        if(ispositive)
            for j=5^str:-1:1
                [tlb,twit] = sk_iterate(X,k,dim,max(tol^2,eps^(3/4)));
                if(tlb >= lb)
                    lb = tlb;
                    lwit = twit;

                    % break out of the function if the target value has already been met
                    if(MetTarget(lb,ub,op_norm,tol,target))
                        lb = op_norm*lb;
                        ub = op_norm*ub;
                        return
                    end
                end
            end
        end
    end

    % Start the semidefinite programming approach for getting upper bounds.
    if(str >= 1 && ((lb + tol < ub && ispositive && nargout > 2) || (ispositive && is_trans_exact && k == 1)))
        cvx_begin sdp quiet
            variable rho(prod_dim,prod_dim) hermitian
            dual variable uwit
            maximize trace(X*rho)
            subject to
                rho >= 0;
                trace(rho) <= 1;
                if(k == 1)
                    uwit : PartialTranspose(rho,2,dim) >= 0;
                else
                    uwit : k*kron(PartialTrace(rho,2,dim),eye(dim(2))) >= rho;
                end
        cvx_end

        % format the output of cvx into a slightly more user-friendly form
        if(strcmpi(cvx_status,'Solved'))
            ub = min(ub,real(cvx_optval));
            if(k == 1)
                uwit = op_norm*PartialTranspose(uwit,2,dim);
            else
                uwit = op_norm*(k*kron(PartialTrace(uwit,2,dim),eye(dim(2))) - uwit);
            end
        else
            error('SkOperatorNorm:NumericalProblems',strcat('Numerical problems encountered (cvx: ',cvx_status,').'));
        end
        
        % In small dimensions, the transpose map gets the result exactly.
        if(is_trans_exact && k == 1 && str >= 1)
            lb = ub;
            [lwit,twit] = eig(rho);
            [~,twit] = max(diag(twit));
            lwit = lwit(:,twit(1));
        elseif(k == 1)
             % we can also get decent lower bounds from the SDP results when k=1
            gs = min(1-roots(jacobi_poly(dim(2)-2,1,1)));
            xmineig = min(real(eig(X)));
            tlb = real(cvx_optval)*(1 - dim(2)*gs/(2*dim(2)-1)) + xmineig*gs/(2*dim(2)-2);
            
            if(tlb > lb)
                lb = tlb;
                lwit = 0; % unfortunately, we don't have a lower bound witness anymore
            end

            % done the str = 1 SDP, now get better upper bounds via symmetric
            % extensions if str >= 2
            for j = 2:str
                 % break out of the function if the target value has already been met
                if(MetTarget(lb,ub,op_norm,tol,target))
                    lb = op_norm*lb;
                    ub = op_norm*ub;
                    return
                end
                
                sym_dim = [dim(1),dim(2)*ones(1,j)];
                prod_sym_dim = dim(1)*dim(2)^j;
                P = kron(eye(dim(1)),SymmetricProjection(dim(2),j));
                    
                cvx_begin sdp quiet
                    variable rho(prod_sym_dim,prod_sym_dim) hermitian
                    dual variable suwit
                    maximize trace(X*PartialTrace(rho,3:j+1,sym_dim))
                    subject to
                        rho >= 0;
                        trace(rho) <= 1;
                        P*rho*P' == rho;
                        suwit : PartialTranspose(rho,1:ceil(j/2)+1,sym_dim) >= 0;
                cvx_end

                % format the output of cvx into a slightly more user-friendly form
                if(strcmpi(cvx_status,'Solved'))
                    if(real(cvx_optval) < ub)
                        ub = real(cvx_optval);
                        uwit = op_norm*PartialTranspose(suwit,1:ceil(j/2)+1,sym_dim);
                    end
                    
                    gs = min(1-roots(jacobi_poly(dim(2)-2,mod(j,2),floor(j/2)+1)));
                    tlb = real(cvx_optval)*(1 - dim(2)*gs/(2*dim(2)-1)) + xmineig*gs/(2*dim(2)-2);

                    if(tlb > lb)
                        lb = tlb;
                        lwit = 0; % unfortunately, we don't have a lower bound witness anymore
                    end
                else
                    error('SkOperatorNorm:NumericalProblems',strcat('Numerical problems encountered (cvx: ',cvx_status,').'));
                end
            end
        end
    end

    lb = op_norm*lb;
    ub = op_norm*ub;
end


% This function checks whether or not the lower bound or upper bound
% already computed meets the desired target value (within numerical error)
% and thus we can abort the script early.
function mt = MetTarget(lb,ub,op_norm,tol,target)
    mt = (op_norm*(lb + tol) >= op_norm*ub || (target >= 0 && (op_norm*(lb - tol) >= target || op_norm*(ub + tol) <= target)));    
end