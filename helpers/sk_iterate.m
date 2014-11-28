%%  SK_ITERATE    Computes a lower bound of the S(k)-norm of an operator
%   This function has one required argument:
%     X: a square positive semidefinite matrix
%
%   SK = sk_iterate(X) is a lower bound of the S(1)-norm of X (i.e., the
%   maximum overlap of X with a separable pure state -- see references
%   [1,2,3]). The two subsystems on which X acts are assumed to have equal
%   dimension in this case (specify the optional DIM parameter if they are
%   of unequal dimension). The lower bound is computed via a randomized
%   method that searches for local maxima. X can be full or sparse.
%
%   This function has four optional arguments:
%     K (the Schmidt rank to optimize, default 1)
%     DIM (default has both subsystems of equal dimension)
%     TOL (default 10^-5)
%     V0 (default is a randomly-chosen vector)
%   
%   [SK,V] = sk_iterate(X,K,DIM,TOL,V0) computes a lower bound of the
%   S(K)-norm of X, as above. The search for a local maximum starts with
%   the vector V0, and the dimensions of the subsystems that X acts on are
%   provided by the 1-by-2 vector DIM. The algorithm terminates when two
%   iterations result in values that are within TOL of each other. V is the
%   local maximizing vector of Schmidt rank <= K.
%
%   URL: http://www.qetlab.com/sk_iterate
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

%   requires: iden.m, MaxEntangled.m, normalize_cols.m, opt_args.m,
%             PermuteSystems.m, SchmidtDecomposition.m, SchmidtRank.m,
%             sporth.m, Swap.m
%             
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: November 14, 2014

%% WHEN UPDATING FOR NON-HERMITIAN, ALSO UPDATE HOW IT'S USED IN SKOPERATORNORM

function [Sk,v] = sk_iterate(X,varargin)

dX = length(X);
sdX = round(sqrt(dX));

% Set optional argument defaults: k=1, dim=sqrt(length(X)), tol=10^-5, v0
% is a random initial vector (set to -1 for now).
[k,dim,tol,v0] = opt_args({ 1, sdX, 10^(-5), -1 },varargin{:});

% allow the user to enter a single number for dim
if(length(dim) == 1)
    dim = [dim,dX/dim];
    if abs(dim(2) - round(dim(2))) >= 2*dX*eps
        error('sk_iterate:InvalidDim','If DIM is a scalar, X must be square and DIM must evenly divide length(X); please provide the DIM array containing the dimensions of the subsystems.');
    end
    dim(2) = round(dim(2));
end
da = dim(1);
db = dim(2);

% Some of this preparation is unnecessary when k = 1, but it's cheap so we
% don't really care.
psi = MaxEntangled(k,1,0);
psiI = Swap(kron(psi,speye(da*db)),[2,3],[k,k,da,db],1);
PSI = psiI*psiI';
PSIX = Swap(kron(psi*psi',X),[2,3],[k,k,da,db],0);

% If the user specified a starting guess v0, parse it; otherwise randomly
% generate one.
randv0 = 1;
if(max(size(v0)) > 1)
    [s0,a0,b0] = SchmidtDecomposition(v0,dim);
    sr = length(s0);
    if(sr > k)
        warning('sk_iterate:SchmidtRankMismatch','The Schmidt rank of the initial vector v0 is %d, which is larger than k=%d. Using a randomly-generated intial vector instead.',sr,k);
    else
        randv0 = 0;
        vp(:,1) = padarray(reshape(a0*diag(s0),da*sr,1),da*(k-sr),'post');
        vp(:,2) = padarray(reshape(b0,db*sr,1),db*(k-sr),'post');
    end
end
if randv0 % generate a random starting vector v0, if appropriate
    vp = randn(max(dim)*k,2) + 1i*randn(max(dim)*k,2);
    v = psiI'*kron(vp(1:k*da,1),vp(1:k*db,2));
    v = v/norm(v);
else
    v = v0;
end
vp = normalize_cols(vp);

% Preparation is done; now do the actual iteration.
it_err = 2*tol+1;
Sk = v'*X*v;
while it_err > tol
    it_err = 0;

    % Loop through the 2 parties.
    for p = 0:1
        % If Schmidt rank is not full, we will have numerical problems; go to
        % lower Schmidt rank iteration.
        sr = SchmidtRank(v,dim,1000*eps);
        if(sr < k)
            [Sk,v] = sk_iterate(X,sr,dim,tol,v);
            break;
        end

        % Fix one of the parties and optimize over the other party.
        if(p==1)
            V0 = kron(vp(1:da*k,1),speye(db*k));
        else
            V0 = kron(speye(da*k),vp(1:db*k,2));            
        end
        V1 = V0'*PSI*V0;

        try
            [vp(1:size(V0,2),p+1),NSk] = eigs(V0'*PSIX*V0,V1,1,'LR');
        catch err
            % In case of ARPACK errors, try to compute the maximal
            % eigenvalue by converting to full and using pinv, if it is
            % reasonable to do so. Otherwise, just ignore it altogether and
            % return the best lower bound found so far.
            if(strcmpi(err.identifier,'MATLAB:eigs:ARPACKroutineError'))
                %if(sdX <= 50)
                %    [tmp,NSk] = eig(pinv(full(V1))*full(V0'*PSIX*V0));
                %    [NSk,ind] = max(diag(NSk));
                %    vp(:,p+1) = tmp(:,ind);
                %else
                    return
                %end
            else
                rethrow(err);
            end
        end
        vp(:,p+1) = vp(:,p+1)/norm(vp(:,p+1));

        it_err = it_err + abs(Sk-NSk);
        Sk = real(NSk);
        v = psiI'*kron(vp(1:k*da,1),vp(1:k*db,2));
        v = v/norm(v);
    end
end