%%  ISBLOCKPOSITIVE    Determines whether or not an operator is block positive
%   This function has one required input argument:
%     X: a square matrix
%
%   IBP = IsBlockPositive(X) is either -1, 0, or 1. A value of 1 indicates
%   that X is block positive, a value of 0 indicates that X is not block
%   positive, and a value of -1 indicates that the block positivity of X
%   could not be determined.
%
%   This function has four optional input arguments:
%     K (default 1)
%     DIM (default has both subsystems of equal dimension)
%     STR (default 2)
%     TOL (default eps^(3/8))
%
%   [IBP,WIT] = IsBlockPositive(X,K,DIM,STR,TOL) is as above, where DIM is
%   a 1-by-2 vector containing the dimensions of the subsystems on which X
%   acts.
%
%   K is the "level" of positivity -- the script checks whether or not X is
%   K-block positive. That is, it checks whether or not X has positive
%   expectation with all states of Schmidt rank <= K.
%
%   STR is an integer that determines how hard the script should work to
%   determine block positivity before giving up (STR = -1 means that the 
%   script won't stop working until it finds an answer). Other valid values 
%   are 0, 1, 2, 3, ... In practice, if STR >= 4 then most computers will
%   run out of memory and/or the sun will explode before computation
%   completes.
%
%   TOL is the numerical tolerance used throughout the script.
%
%   WIT is a witness that verifies that X is (or is not) block positive.
%
%   URL: http://www.qetlab.com/IsBlockPositive

%   requires: cvx (http://cvxr.com/cvx/), iden.m, IsPSD.m, kpNorm.m,
%             MaxEntangled.m, normalize_cols.m, opt_args.m, PartialMap.m,
%             PartialTrace.m, PartialTranspose.m, PermuteSystems.m,
%             Realignment.m, SchmidtDecomposition.m, SchmidtRank.m,
%             sk_iterate.m, SkOperatorNorm.m, SkVectorNorm.m, sporth.m,
%             Swap.m
%
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: September 23, 2014

function [ibp,wit] = IsBlockPositive(X,varargin)

    wit = 0;
    X = full(X);
    dX = size(X);

    % set optional argument defaults: k=1, dim=sqrt(length(X)), str=2,
    % tol=eps^(3/8)
    [k,dim,str,tol] = opt_args({ 1, round(sqrt(dX(1))), 2, eps^(3/8) },varargin{:});
    if(str == -1)
        str = 1/eps; % keep going forever!
    end

    % Make sure that X is Hermitian.
    if(dX(1) ~= dX(2) || (max(max(abs(X-X'))) > tol))
        ibp = 0;
        return
    end
    
    % allow the user to enter a single number for dim
    if(length(dim) == 1)
        dim = [dim,dX(1)/dim];
        if abs(dim(2) - round(dim(2))) >= 2*dX(1)*eps
            error('IsBlockPositive:InvalidDim','If DIM is a scalar, it must evenly divide length(X); please provide the DIM array containing the dimensions of the subsystems.');
        end
        dim(2) = round(dim(2));
    end
    
    % When a local dimension is small, block positivity is trivial.
    if(min(dim) <= k)
        ibp = IsPSD(X);
        wit = zeros(dX);
        return        
    end
    
    op_norm = norm(X);
    Y = op_norm*speye(dX(1)) - X; % we compute the S(k)-norm of this operator
    
    [lb,lwit,ub,uwit] = SkOperatorNorm(Y,k,dim,str,op_norm,tol); % compute the norm

    if(ub <= op_norm*(1 + tol)) % block positive
        ibp = 1;
        wit = uwit;
    elseif(lb >= op_norm*(1 - tol)) % not block positive
        ibp = 0;
        wit = lwit;
    else % not sure :(
        ibp = -1;
    end
end