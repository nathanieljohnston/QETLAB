%%  SYMMETRICEXTENSION    Determines whether or not an operator has a symmetric extension
%   This function has one required argument:
%     X: a positive semidefinite matrix
%
%   EX = SymmetricExtension(X) is either 1 or 0, indicating that X does or
%   does not have a 2-copy symmetric extension. The extension is always
%   taken on the second subsystem of X.
%
%   This function has five optional arguments:
%     K (default 2)
%     DIM (default has both subsystems of equal dimension)
%     PPT (default 0)
%     BOS (default 0)
%     TOL (default eps^(1/4))
%
%   [EX,WIT] = SymmetricExtension(X,K,DIM,PPT,BOS,TOL) determines whether
%   or not X has a symmetric extension and provides a witness WIT that
%   verifies the answer. If a symmetric extension of X exists
%   (i.e., EX = 1) then WIT is such a symmetric extension. If no symmetric
%   extension exists (i.e., EX = 0) then WIT is an entanglement witness
%   with trace(WIT*X) = -1 but trace(WIT*Y) >= 0 for all symmetrically
%   extendable Y.
%
%   K is the desired number of copies of the second subsystem. DIM is a
%   1-by-2 vector containing the dimensions of the subsystems on which X
%   acts. PPT is a flag (either 0 or 1) indicating whether the desired
%   symmetric extension must have positive partial transpose. BOS is a flag
%   (either 0 or 1) indicating whether the desired symmetric extension must
%   be Bosonic (i.e., be supported on the symmetric subspace). TOL is the
%   numerical tolerance used when determining whether or not a symmetric
%   extension exists.
%
%   URL: http://www.qetlab.com/SymmetricExtension

%   requires: cvx (http://cvxr.com/cvx/), IsPPT.m, IsPSD.m, opt_args.m,
%             PartialTrace.m, PartialTranspose.m, PermutationOperator.m,
%             PermuteSystems.m, sporth.m, SymmetricProjection.m
%
%   author: Nathaniel Johnston (nathaniel@njohnston.ca), with improvements
%           by Mateus Ara�jo
%   package: QETLAB
%   last updated: January 5, 2023

function [ex,wit] = SymmetricExtension(X,varargin)

lX = length(X);

% set optional argument defaults: k=2, dim=sqrt(length(X)), ppt=0, bos=0, tol=eps^(1/4)
[k,dim,ppt,bos,tol] = opt_args({ 2, round(sqrt(lX)), 0, 0, eps^(1/4) },varargin{:});

% allow the user to enter a single number for dim
if(length(dim) == 1)
    dim = [dim,lX/dim];
    if abs(dim(2) - round(dim(2))) >= 2*lX*eps
        error('SymmetricExtension:InvalidDim','If DIM is a scalar, it must evenly divide length(X); please provide the DIM array containing the dimensions of the subsystems.');
    end
    dim(2) = round(dim(2));
end

% In certain situations, we don't need semidefinite programming.
if(k == 1 || (lX <= 6 && ppt && nargout <= 1))
    if(k == 1 && ~ppt) % in some cases, the problem is *really* trivial
        if(nargout > 1)
            [ex,wit] = IsPSD(X,tol);
            if(ex == 1)
                wit = X;
            end
        else
            ex = IsPSD(X,tol);
        end
        
        return
    end
    
    % In this case, all they asked for is a 1-copy PPT symmetric extension
    % (i.e., they're asking if the state is PPT).
    if(nargout > 1)
        [ex,wit] = IsPPT(X,2,dim,tol);
        if(ex)
            wit = X;
        else
            wit = PartialTranspose(wit*wit');
            wit = -wit/trace(wit*X);
        end
    else
        ex = (IsPPT(X,2,dim,tol) && IsPSD(X,tol));
    end

% In the 2-qubit case, an analytic formula is known for whether or not a
% state has a (2-copy, non-PPT) symmetric extension that is much
% faster to use than semidefinite programming.
elseif(~isa(X,'cvx') && k == 2 && ~ppt && dim(1) == 2 && dim(2) == 2 && nargout <= 1) % we don't need "&& ~bos" thanks to a lemma of Myhr and Lutkenhaus
    ex = (trace(PartialTrace(X,1)^2) >= trace(X^2) - 4*sqrt(det(X)));
    
% otherwise, use semidefinite programming to find a symmetric extension
elseif(k > 1)
    sdp_dim = [dim(1),dim(2)*ones(1,k)];
    % For Bosonic extensions, it suffices to optimize over the symmetric
    % subspace, which is smaller.
    if(bos)
        sdp_prod_dim = dim(1)*nchoosek(k+dim(2)-1, dim(2)-1);
    else
        sdp_prod_dim = dim(1)*dim(2)^k;
    end
    cvx_begin sdp quiet
        cvx_precision(tol);
        if(bos)
        	V = kron(speye(dim(1)),SymmetricProjection(dim(2),k,1)); %this is an isometry
        	variable sig(sdp_prod_dim,sdp_prod_dim) hermitian
        	rho = V*sig*V';
        else
        	variable rho(sdp_prod_dim,sdp_prod_dim) hermitian
        end
        if(nargout > 1) % don't want to compute the dual solution in general (especially not if this is called within CVX)
            dual variable W
        end
        subject to
            if(nargout > 1)
                W : PartialTrace(rho,3:k+1,sdp_dim) == X;
            else
                PartialTrace(rho,3:k+1,sdp_dim) == X;
            end
            if(ppt)
            	for j=2:k+1
                	PartialTranspose(rho,2:j,sdp_dim) >= 0;
                end
            end
        	if(bos)
        		sig >= 0;
        	else
        		rho >= 0;
                for j = 3:k+1% in the permutation invariant case we need to explicitly enforce the symmetry
                    PartialTrace(rho,setdiff(2:k+1,j),sdp_dim) == X;
                end
        	end
    cvx_end
    
    ex = 1-min(cvx_optval,1); % CVX-safe way to map (0,Inf) to (1,0)
    if(~isa(X,'cvx')) % make the output prettier if it's not a CVX input
        ex = round(ex);

        % Deal with error messages and witnesses.
        if(nargout > 1 && strcmpi(cvx_status,'Solved'))
            wit = rho;
        elseif(strcmpi(cvx_status,'Inaccurate/Solved'))
            wit = rho;
            warning('SymmetricExtension:NumericalProblems','Minor numerical problems encountered by CVX. Consider adjusting the tolerance level TOL and re-running the script.');
        elseif(nargout > 1 && strcmpi(cvx_status,'Infeasible'))
            wit = -W;
        elseif(strcmpi(cvx_status,'Inaccurate/Infeasible'))
            if(nargout > 1)
                wit = -W;
            end
            warning('SymmetricExtension:NumericalProblems','Minor numerical problems encountered by CVX. Consider adjusting the tolerance level TOL and re-running the script.');
        elseif(strcmpi(cvx_status,'Unbounded') || strcmpi(cvx_status,'Inaccurate/Unbounded') || strcmpi(cvx_status,'Failed'))
            error('SymmetricExtension:NumericalProblems',strcat('Numerical problems encountered (CVX status: ',cvx_status,'). Please try adjusting the tolerance level TOL.'));
        end
    end
else
    error('SymmetricExtension:InvalidK','K must be a positive integer.');
end
