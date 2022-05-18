%%  RANDOMPPTSTATE    Generates a random density matrix with positive partial transpose, and optionally low rank
%   This function has one required argument:
%     DIM: a scalar or 2-entry vector specifying the local dimenions of the
%          desired PPT density matrix
%
%   RHO = RandomPPTState(DIM) generates a density matrix with local
%   dimensions DIM(1) and DIM(2) (or both local dimensions equal to DIM, if
%   DIM is a scalar) with the property that its partial transpose is
%   positive semidefinite. The random generation does not follow any
%   particular well-known or named distribution.
%
%   This function has three optional arguments:
%     RNK (default DIM(1)*DIM(2))
%     TOL (default 10^(-12))
%     MAX_ITS (default infinity)
%
%   RHO = RandomPPTState(DIM,RNK,TOL,MAX_ITS) generates a random PPT
%   density matrix with the property that its rank is at most RNK(1) and
%   the rank of its partial transpose is at most RNK(2) (RNK can optionally
%   be provided as a scalar, in which case both ranks are bounded by RNK).
%   This state is computed via the algorithm of "Low-rank extremal
%   positive-partial-transpose states and unextendible product bases" by
%   Leinass et al. TOL is a numerical tolerance used in the calculation
%   (all eigenvalues that we wish were zero are guaranteed to have absolute
%   value smaller than TOL). MAX_ITS is the maximum number of iterations to
%   carry out; set this to some finite value if you want long-running
%   calculations to end early.
%
%   URL: http://www.qetlab.com/RandomPPTState

%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: May 18, 2022

% REF: https://arxiv.org/pdf/1002.1949.pdf

function rho = RandomPPTState(dim,varargin)
    % allow the user to enter a single number for dim
    if(length(dim) == 1)
        dim = [dim,dim];
    end
    pDim = prod(dim);
    
    % set optional argument defaults: rnk=[pDim,pDim],tol=10^(-12),max_its=Inf
    [rnk,tol,max_its] = opt_args({ pDim, 10^(-12), Inf },varargin{:});
    
    % allow the user to enter a single number for rnk
    if(length(rnk) == 1)
        rnk = [rnk,rnk];
    end

    % For implementation reasons, we want either both ranks to be full, or
    % neither; not one or the other.
    if(rnk(1) ~= rnk(2) && rnk(1) >= pDim)
        rnk(1) = min(rnk(1),pDim-1);
    elseif(rnk(1) ~= rnk(2) && rnk(2) >= pDim)
        rnk(2) = min(rnk(2),pDim-1);
    end

    % Generate a random starting state.
    rho = RandomDensityMatrix(pDim);

    % Shift it so that it has PPT.
    lam_min = abs(min(real(eig(PartialTranspose(rho,2,dim)))));
    rho = rho + lam_min*eye(pDim);
    rho = rho/trace(rho);
   
    % If both ranks can be full, just be done.
    if(min(rnk) >= pDim)
        return;
    end
    
    % Otherwise, iterate!
    ct = 0;
    while 1
        [rho,kminev] = PPTIterate(rho,dim,rnk);
        rho = rho + rho';
        rho = rho/trace(rho);
        
        % If the desired tolerance has been reached, or the maximum number of iterations has been reached, be done.
        ct = ct + 1;
        if(kminev < tol || ct >= max_its)
            break;
        end
    end
end

function [sigma,kminev] = PPTIterate(rho,dim,k)
    pDim = prod(dim);
    [vrho,drho] = eig(rho,'vector');
    [vrhoPT,drhoPT] = eig(PartialTranspose(rho,2,dim),'vector');
    
    % Get rid of potentialy tiny imaginary pieces of eigenvalues, resulting
    % from numerical precision issues.
    drho = real(drho);
    drhoPT = real(drhoPT);
    
    % Make sure eigenvalues are sorted, and keep eigenvectors sorted in the
    % same way.
    [drho,ind] = sort(drho,'descend');
    vrho = vrho(:,ind);
    [drhoPT,ind] = sort(drhoPT,'descend');
    vrhoPT = vrhoPT(:,ind);
    
    kminev = max(abs(drho(k(1)+1)),abs(drhoPT(k(2)+1)));
    
    % Construct the matrix B and vector mu used in the iteration.
    mu = [drho(k(1)+1:end);drhoPT(k(2)+1:end)];
    lam = [vrho(:,k(1)+1:end),vrhoPT(:,k(2)+1:end)]; % matrix whose columns are eigenvectors associated with eigenvalues in mu

    B = zeros(2*pDim-k(1)-k(2),pDim^2);
    for j = 1:(pDim-k(1))
        B(j,:) = sym_vectorize(lam(:,j)*lam(:,j)')';
    end
    for j = (pDim-k(1)+1):(2*pDim-k(1)-k(2))
        B(j,:) = sym_vectorize(PartialTranspose(lam(:,j)*lam(:,j)',2,dim))';
    end
    
    del_x = -pinv(B)*mu;
    new_x = sym_vectorize(rho) + del_x;
    sigma = unsym_vectorize(new_x);
end

% Compute an orthonormally-scaled symmetric vectorization of a Hermitian
% matrix.
function sv = sym_vectorize(X)
    n = length(X);
    sX = diag(diag(X)) + sqrt(2)*(real(triu(X,1)) + imag(tril(X,-1)));
    sv = reshape(sX,n^2,1);
end

% Inverse function of sym_vectorize.
function X = unsym_vectorize(symX)
    n = round(sqrt(length(symX)));
    tmpX = reshape(symX,n,n);
    triX = (triu(tmpX,1) + 1i*tril(tmpX,-1))/sqrt(2);
    X = diag(diag(tmpX)) + triX + triX';
end