%%  DISTINGUISHABILITY    Computes the maximum probability of distinguishing quantum states
%   This function has one required input argument:
%     X: a cell containing density matrices or a single matrix containing
%        pure states as its column vectors
%
%   DIST = Distinguishability(X) is the maximum probability of
%   distinguishing the quantum states specified by X. X can either be a
%   cell containing 2 or more density matrices to be distinguished, or X
%   can be a matrix whose columns are pure vector states to be
%   distinguished.
%
%   This function has one optional input argument:
%     P (default [1/k, ..., 1/k], where k is the number of states)
%
%   [DIST,MEAS] = Distinguishability(X,P) is the maximum probability of
%   distinguishing the quantum states specified by X, where the vector P
%   contains the probability that each state is chosen (by default, the
%   states are chosen uniformly at random). MEAS is a cell containing the
%   optimal measurement operators.
%
%   URL: http://www.qetlab.com/Distinguishability

%   requires: cvx (http://cvxr.com/cvx/), kpNorm.m, opt_args.m, TraceNorm.m
%
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: October 6, 2014

function [dist,meas] = Distinguishability(X,varargin)

if(iscell(X))
    num_ops = length(X);
    dim = length(X{1});
else
    [dim,num_ops] = size(X);
end
    
% set optional argument defaults: p = [1/k,1/k,...,1/k], where k = # of states
[p] = opt_args({ ones(1,num_ops)/num_ops },varargin{:});

if(abs(sum(p) - 1) > num_ops^2*eps || length(p) ~= num_ops)
    error('Distinguishability:InvalidP','The vector P must be a probability distribution of the same length as the number of states: its elements must be non-negative and they must sum to 1.');
end

if(num_ops == 1 || max(p) >= 1) % of course we can distinguish 1 object
    dist = 1;
    if(nargout > 1) 
        meas = eye(dim); % optimal measurements is trivial in this case
    end
    return
end

% X is a cell of density matrices
if(iscell(X))
    for j = 1:num_ops % make sure that the density operators are scaled
        X{j} = X{j}/trace(X{j});
    end
    
    % There is a closed-form expression for the distinguishability of two
    % density matrices.
    if(num_ops == 2)
        dist = 1/2 + TraceNorm(p(1)*X{1} - p(2)*X{2})/2;
        if(nargout > 1) % construct optimal measurements, if requested
            [v,d] = eig(p(1)*X{1} - p(2)*X{2});
            d = diag(d);
            pind = find(d >= 0);
            meas{1} = v(:,pind)*v(:,pind)';
            meas{2} = eye(dim) - meas{1};
        end
        return
    end
    
    % Check to see if the states are mutually orthogonal, and return 1 if
    % they are.
    if(num_ops <= dim)
        mut_orth = 1;
        for j = 1:num_ops
            for k = j+1:num_ops
                if(max(max(abs(X{j}*X{k}))) > eps*dim^2)
                    mut_orth = 0;
                    break;
                end
            end
        end
        if(mut_orth == 1) % states are mutually orthogonal
            dist = 1;
            if(nargout > 1) % construct optimal measurements, if requested
                meas_sum = zeros(dim);
                for j = num_ops:-1:2 % pre-allocate for speed
                    oX = orth(X{j});
                    meas{j} = oX*oX';
                    meas_sum = meas_sum + meas{j};
                end
                meas{1} = eye(dim) - meas_sum;
            end
            return;
        end
    end
    
% X is a matrix whose columns are pure states
else
    X = normalize_cols(X); % make sure that the columns have unit length
    
    % There is a closed-form expression for the distinguishability of two
    % pure states.
    if(num_ops == 2)
        dist = 1/2 + sqrt(2*(p(1)^2 + p(2)^2) - 4*p(1)*p(2)*abs(X(:,1)'*X(:,2))^2)/2;
        if(nargout > 1) % construct optimal measurements, if requested
            [v,d] = eig(p(1)*X(:,1)*X(:,1)' - p(2)*X(:,2)*X(:,2)');
            d = diag(d);
            pind = find(d >= 0);
            meas{1} = v(:,pind)*v(:,pind)';
            meas{2} = eye(dim) - meas{1};
        end
        return
    end
    
    % Check to see if the states are mutually orthogonal, and return 1 if
    % they are.
    if(num_ops <= dim)
        X2 = X'*X;
        if(max(max(abs(X2 - diag(diag(X2))))) < eps*dim^2)
            dist = 1; % the states are orthogonal, so perfectly distinguishable
            if(nargout > 1) % construct optimal measurements, if requested
                meas_sum = zeros(dim);
                for j = num_ops:-1:2 % pre-allocate for speed
                    meas{j} = X(:,j)*X(:,j)';
                    meas_sum = meas_sum + meas{j};
                end
                meas{1} = eye(dim) - meas_sum;
            end
            return
        end
    end
    
    % Turn pure states into density operator form.
    Y = X; X = {};
    for j = num_ops:-1:1 % pre-allocate for speed
        X{j} = Y(:,j)*Y(:,j)';
    end
end

% For 3 or more density matrices, we have to use semidefinite programming.
cvx_begin sdp quiet
    cvx_precision default;
    variable P(dim,dim,num_ops) hermitian
    
    P_tr = 0;
    P_sum = zeros(dim);
    for j = 1:num_ops
        P_tr = P_tr + p(j)*trace(P(:,:,j)*X{j});
        P_sum = P_sum + P(:,:,j);
    end
    P_tr = P_tr + P_tr';
    
    maximize P_tr
    subject to
        P_sum == eye(dim);
        P >= 0;
cvx_end
dist = real(cvx_optval)/2;

% Also return the optimal measurements, if requested.
if(nargout > 1)
    meas = mat2cell(reshape(P,dim,dim*num_ops),dim,dim*ones(1,num_ops));
end