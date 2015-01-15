%%  ISABSPPT    Determines whether or not a density matrix is absolutely PPT
%   This function has one required input argument:
%     RHO: a density matrix, or a vector of the eigenvalues of a density
%          matrix
%
%   IAPPT = IsAbsPPT(RHO) is either 1 or 0, indicating that RHO (which is
%   assumed to act on bipartite space with each of the two subsystems of
%   equal dimension) is or is not absolutely PPT (see [1] for the
%   definition of absolutely PPT states). If the two subsystems are both
%   of dimension 7 or higher, than this function may sometimes return a
%   value of -1, indicating that it could not determine whether or not RHO
%   is absolutely PPT within a reasonable amount of time.
%
%   This function has on optional arguments:
%     DIM (default has both subsystems of equal dimension)
%
%   IAPPT = IsAbsPPT(RHO,DIM) is as above, where DIM is a 1-by-2 vector
%   containing the dimensions of the subsystems on which RHO acts.
%
%   URL: http://www.qetlab.com/IsAbsPPT
%
%   References:
%   [1] R. Hildebrand. Positive partial transpose from spectra. Phys. Rev.
%       A, 76:052325, 2007. E-print: arXiv:quant-ph/0502170

%   requires: AbsPPTConstraints.m, InSeparableBall.m, IsPSD.m, opt_args.m,
%             perm_inv.m
%
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: December 24, 2014

function iappt = IsAbsPPT(rho,varargin)

iappt = -1; % by default, we don't know the answer
sz = size(rho);
len = max(sz);

% set optional argument defaults: dim=sqrt(length(rho))
[dim] = opt_args({ round(sqrt(len)) },varargin{:});

% If rho is a vector, it is a vector of eigenvalues. If it is a matrix,
% we should find its eigenvalues.
if(min(sz) == 1) % rho is a vector of eigenvalues
    if(~isa(rho,'cvx'))
        lam = sort(real(rho),'descend');
    else
        lam = rho;
    end
else % rho is a density matrix
    if(isa(rho,'cvx'))
        error('IsAbsPPT:InvalidRHO','If used within a CVX optimization problem, RHO must be a vector of eigenvalues, not a density matrix.');
    else
        lam = sort(real(eig(full(rho))),'descend');
    end
end

% allow the user to enter a single number for dim
if(length(dim) == 1)
    dim = [dim,len/dim];
    if abs(dim(2) - round(dim(2))) >= 2*len*eps
        error('IsAbsPPT:InvalidDim','If DIM is a scalar, it must evenly divide length(RHO); please provide the DIM array containing the dimensions of the subsystems.');
    end
    dim(2) = round(dim(2));
end
pd = prod(dim);
p = min(dim);

% Now do some error checking.
if(pd ~= length(lam))
    error('IsAbsPPT:InvalidDim','The dimensions provided by DIM do not match the length of RHO.');
end

% Do fast checks first, if the input is numeric.
if(~isa(lam,'cvx'))
    % Check if it's in the Gurvits-Barnum ball.
    if(InSeparableBall(lam))
        iappt = 1;
        return;

    % Check if it satisfies the Gerschgorin circle property from arXiv:1406.1277.
    elseif(sum(lam(1:p-1)) <= 2*lam(end) + sum(lam(end-p+1:end-1)))
        iappt = 1;
        return;
    end
end

% Still don't know the answer, and dim is small enough? Roll up your
% sleeves and *really* check whether or not it's absolutely PPT.
if(isa(lam,'cvx'))
    if(p >= 7)
        warning('IsAbsPPT:LargeDim','Only optimizing over the first 2612 linear matrix inequalities for absolutely PPT states, since there are over a million constraints in this case. Thus the set being optimized over is slightly larger than the true set of absolutely PPT states.');
    end
    
    cvx_begin sdp quiet
    L = AbsPPTConstraints(lam,dim,0,2612);
    subject to
        for j = 1:length(L)
            L{j} >= 0;
        end
        for j = 1:pd-1
            lam(j) >= lam(j+1); % eigenvalues must be in order
        end
        lam >= 0;
    cvx_end
    iappt = 1-min(cvx_optval,1); % CVX-safe way to map (0,Inf) to (1,0)
else
    % Compute *all* constraints up to 6 local dimensions, otherwise only
    % compute a few thousand constraints.
    L = AbsPPTConstraints(lam,dim,1,2612);
    if(~IsPSD(L{end})) % only need to check this one constraint: the AbsPPTConstraints function already checked all the earlier ones
        iappt = 0;
    elseif(p <= 6) % the test we just did is complete for local dimensions up to 6
        iappt = 1;
    end
end