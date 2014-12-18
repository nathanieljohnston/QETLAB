function [DIST, MEAS, DUAL_SOL] = LocalDistinguishability(X, varargin)
%%  LocalDistinguishability   Computes the maximum probability of distinguishing quantum states by symmetric-extendible measurements
%   This function has one required argument:
%      X: a cell array containing density matrices, or a single matrix 
%           containing pure states as its column vectors
%   
%   This function has six optional arguments:
%     P (default [1/N, ..., 1/N], where N is the number of states)
%     DIM (default has both subsystems of equal dimension)
%     COPIES (default 2)
%     PPT (default true)
%     BOS (default true)
%     TOL (default eps^(1/4))
%
%   [DIST, MEAS, DUAL_SOL] = LocalDistinguishability(X) returns the optimal 
%      probability to distinguish the input states by a measurement whose 
%      operators have 2-copies PPT bosonic symmetric extensions. MEAS is a 
%       cell array containing the optimal measurement operators. DUAL_SOL
%       is the optimal solution of the dual problem.
%
%   COPIES is the desired number of copies of the second subsystem. 
%   DIM is a 1-by-2 vector containing the dimensions of the subsystems 
%   on which X acts. 
%   PPT is a flag (either true or false) indicating whether the symmetric
%   extensions of the measurement operators must have positive partial 
%   transpose. 
%   BOS is a flag (either true or false) indicating whether the symmetric
%   extensions of the measurement operators must be Bosonic (i.e., be 
%   supported on the symmetric subspace). 
%   TOL is the numerical tolerance used when determining whether or not 
%   a symmetric extension exists.
%
%   The arguments are mutually optional (note this is different from other 
%   QETLAB functions!). For example, if you want to specify the number of 
%   COPIES to 3 and the dimension of the two subsystems to [2 3], you can 
%   call the function as follows:
%   LocalDistinguishability(X, 'COPIES', 3, 'DIM', [2 3])
%
%   URL: http://www.qetlab.com/LocalDistinguishability
%
%   author: Alessandro Cosentino (cosenal@gmail.com)

if(iscell(X))
    num_ops = length(X);
    d = length(X{1});
    S = X;
    for j = 1:num_ops % make sure that the density operators are scaled
        S{j} = S{j}/trace(S{j});
    end
else
    X = normalize_cols(X);  % make sure that the columns have unit length
    [d,num_ops] = size(X);
    S = mat2cell(X, d, ones(1, num_ops));
    S = cellfun(@(x) x*x', S, 'UniformOutput', false);
end

parser = inputParser;
addRequired(parser, 'X', @isnumeric);
addParameter(parser, 'P', ones(1,num_ops)/num_ops, @isnumeric);
addParameter(parser, 'DIM', round(sqrt(d)), @isnumeric);
addParameter(parser, 'COPIES', 2, @isscalar);
addParameter(parser, 'PPT', true, @islogical);    
addParameter(parser, 'BOS', true, @islogical);
addParameter(parser, 'TOL', eps^(1/4), @isscalar);
parse(parser,X,varargin{:});
p = parser.Results.P;
dim = parser.Results.DIM;
copies = parser.Results.COPIES;
ppt = int8(parser.Results.PPT);
bos = int8(parser.Results.BOS);
tol = parser.Results.TOL;

if(abs(sum(p) - 1) > num_ops^2*eps || length(p) ~= num_ops)
    error('LocalDistinguishability:InvalidP', ...
        strcat('The vector P must be a probability distribution of ', ...
        'the same length as the number of states: its elements must', ...
        'be non-negative and they must sum to 1.'));
end

if(num_ops == 1 || max(p) >= 1) % of course we can distinguish 1 object
    DIST = 1;
    if(nargout > 1) 
        MEAS = eye(dim); % optimal measurements is trivial in this case,
        DUAL_SOL = S{1}; % and so is the solution of the dual problem
    end
else
    cvx_begin sdp quiet
        cvx_precision default;
        variable P(d,d,num_ops) hermitian
        dual variable Y

        P_tr = 0;
        P_sum = zeros(d);
        for j = 1:num_ops
            P_tr = P_tr + p(j)*trace(P(:,:,j)*S{j});
            P_sum = P_sum + P(:,:,j);
        end
        P_tr = P_tr + P_tr';

        maximize P_tr
        subject to
            Y : eye(d) == P_sum;
            for j = 1:num_ops
                SymmetricExtension(P(:,:,j),copies,dim,ppt,bos,tol) == 1;
            end
    cvx_end

    DIST = real(cvx_optval)/2;
    
    % Also return the optimal measurements and the optimal solution
    % of the dual problem, if requested.
    if(nargout > 1)
        MEAS = mat2cell(reshape(P,d,d*num_ops),d,d*ones(1,num_ops));
        DUAL_SOL = Y;
    end
end
end