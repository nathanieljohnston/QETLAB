function [DIST, MEAS, DUAL_SOLUTION] = LocalDistinguishability(X, varargin)
%%  LocalDistinguishability
%   This function has one required argument:
%      X: a cell array containing density matrices, or a single matrix 
%           containing pure states as its column vectors
%
%   v = LocalDistinguishability(X, 'copies', c) returns the optimal 
%     probability to distinguish the input states by a measurement whose 
%     operators have c-copies symmetric extensions.
%   
%   This function has four optional arguments:
%     p (default [1/N, ..., 1/N], where N is the number of states)
%     dim (default has both subsystems of equal dimension)
%     copies (default 1)
%     ppt (default true)
%
%   The arguments are mutually optional. For example, if you want to
%   specify the number of copies to 2 and the dimension of the two 
%   subsystems to [2 3], you call the function with the following 
%   arguments:
%   LocalDistinguishability(X, 'copies', 2, 'dim', [2 3])
%
%   URL: http://www.qetlab.com/LocalDistinguishability
%
%   requires: cvx (http://cvxr.com/cvx/), PartialTrace.m, 
%             PartialTranspose.m, PermutationOperator.m,
%             PermuteSystems.m, sporth.m, SymmetricProjection.m
%
%   author: Alessandro Cosentino (cosenal@gmail.com)

    if(iscell(X))
        N = length(X);
        d = length(X{1});
        S = X;
    else
        [d,N] = size(X);
        S = mat2cell(X, d, ones(1, N));
        S = cellfun(@(x) x*x', S, 'UniformOutput', false);
    end
    
    parser = inputParser;
    addRequired(parser, 'X', @isnumeric);
    addParameter(parser, 'p', ones(1,N)/N, @isnumeric);
    addParameter(parser, 'dim', round(sqrt(d)), @isnumeric);
    addParameter(parser, 'copies', 1, @isscalar);
    addParameter(parser, 'ppt', true, @islogical);    
    parse(parser,X,varargin{:});
    p = parser.Results.p;
    dim = parser.Results.dim;
    copies = parser.Results.copies;
    ppt = parser.Results.ppt;
    
    if(isscalar(dim))
        dim = [dim,d/dim];
        if (abs(dim(2) - round(dim(2))) >= 2*d*eps)
            error('SymmetricExtDist:InvalidDim', 'If DIM is a scalar, it must evenly divide length(X); please provide the dim array containing the dimensions of the subsystems.');
        end
        dim(2) = round(dim(2));
    end
    
    if (copies > 3)
        warning('SymmetricExtensionDist:ManyCopies', ...
            'Checking %d-copies symmetric extensions takes time!', copies);
    end
        
    if(length(p) ~= N || sum(p) ~= 1)
         error('SymmetricExtensionDist:InvalidProb', ...
             'Invalid probability vector');
    else
        for k=1:N
            S{k} = S{k}*p(k);    
        end
    end
        
    Proj = kron(speye(dim(1)), SymmetricProjection(dim(2), copies, 1));
    sdp_prod_dim = size(Proj, 2);
    sdp_dim = [dim(1), dim(2)*ones(1,copies)];
    
    cvx_begin sdp quiet
        cvx_precision default

        variable X(sdp_prod_dim, sdp_prod_dim, N) hermitian
        dual variable Y

        opt = 0;
        P_sum = zeros(d);
        P = cell(1, N);
        for k=1:N
            P{k} = PartialTrace(Proj*X(:,:,k)*Proj', 3:copies+1, sdp_dim);
            opt = opt + trace(P{k}*S{k});
            P_sum = P_sum + P{k};
        end 

        maximize(opt);
        subject to
            Y : eye(d) == P_sum;
            for k=1:N
                if (ppt)
                    PartialTranspose(Proj*X(:,:,k)*Proj', ...
                        1:floor(copies/2)+1, sdp_dim) >= 0;
                end
                X(:,:,k) >= 0;
            end        
    cvx_end

    DIST = real(cvx_optval);
    MEAS = P;
    DUAL_SOLUTION = Y;
end
