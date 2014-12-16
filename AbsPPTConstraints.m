%%  ABSPPTCONSTRAINTS    Builds the eigenvalue matrices that determine whether or not a state is absolutely PPT
%   This function has one required input argument:
%     LAM: a vector of length prod(DIM) that contains the eigenvalues of
%          the state
%
%   L = AbsPPTConstraints(LAM) is a cell containing the constraint
%   matrices constructed in [1] that characterize whether or not a mixed
%   state with eigenvalues LAM is "absolutely PPT" (see [1] for
%   definition). In particular, a state with eigenvalues LAM is absolutely
%   PPT if and only if every matrix in L is positive semidefinite.
%
%   These matrices can be used as constraints in a CVX optimization
%   problem, allowing the user to solve SDPs whose feasible region is the
%   set of absolutely PPT states (see online documentation for examples).
%
%   This function has two optional input arguments:
%     DIM (default has both subsystems of equal dimension)
%     ESC_IF_NPOS (default 0)
%     LIM (default 0)
%
%   L = AbsPPTConstraints(LAM,DIM,ESC_IF_NPOS,LIM) is as above, but if LIM
%   > 0 then at most LIM matrices are computed. The two local dimensions
%   are specified by the 1-by-2 vector DIM. Additionally, if ESC_IF_NPOS =
%   1 then the function stops computing as soon as a matrix in L is not
%   positive semidefinite.
%
%   URL: http://www.qetlab.com/AbsPPTConstraints
%
%   References:
%   [1] R. Hildebrand. Positive partial transpose from spectra. Phys. Rev.
%       A, 76:052325, 2007. E-print: arXiv:quant-ph/0502170

%   requires: IsPSD.m, opt_args.m, perm_inv.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: December 2, 2014

function L = AbsPPTConstraints(lam,varargin)

    sz = size(lam);
    if(min(sz) > 1) % LAM is a density matrix: get its eigenvalues
        lam = eig(lam);
    end
    lam_len = max(sz);
    dim_guess = round(sqrt(lam_len));
    
    % set optional argument defaults: dim=dim_guess, esc_if_npos=0, lim=0
    [dim,esc_if_npos,lim] = opt_args({ dim_guess, 0, 0 },varargin{:});

    % Allow the user to enter a single number for DIM.
    if(length(dim) == 1)
        dim = [dim,dim];
    end
    
    % Now do some error checking.
    if(prod(dim) ~= lam_len)
        error('AbsPPTConstraints:InvalidDim','The dimensions provided by DIM do not match the length of the vector LAM.');
    end
    
    % Make sure that LAM is sorted in descending order (if possible).
    if(isa(lam,'cvx') == 0)
        lam = sort(lam,'descend');
    end
    
    % No simple formula is known for the number of matrices in L that will
    % be returned, so just increase L's size every step. The speed
    % reduction from re-allocating memory is minimal compared to the time
    % it takes to do the recursion anyway.
    L = cell(0);
    
    % Begin by setting some preliminary values that will let the search
    % work.
    esc_now = 0;
    p = min(dim);
    X = (p*(p+1)/2+1)*ones(p); % set each entry to anything larger than p*(p+1)/2
    num_pool = ones(p*(p+1)/2,1);
    ct = 1;

    % If p <= 2, hard-code the results (the recursion will get cranky).
    if(p == 1)
        return;
    elseif(p == 2)
        L{1} = [2*lam(end), lam(end-1)-lam(1); lam(end-1)-lam(1), 2*lam(end-2)];
        return;
    end
    
    % The two top-left and bottom-right entries of the matrix are always
    % the same, so fix them.
    X(1,1) = 1;
    X(1,2) = 2;
    X(p,p) = p*(p+1)/2;
    X(p-1,p) = p*(p+1)/2 - 1;

    % Now recursely construct the absolutely PPT constraints.
    fill_matrix(1,3,3);

    % This function does the recursion itself. It is nested because it
    % needs to have access to variables like num_pool and X.
    function fill_matrix(i,j,low_lim)
        up_lim = min(j*(j+1)/2 + i*(p-j), p*(p+1)/2 - 2);
        for k = low_lim:up_lim
            if(num_pool(k) == 1)
                X(i,j) = k;
                num_pool(k) = 0;

                if(check_ordered(X,i,j) == 1)
                    if(i == p-1 && j == p-1) % we have completely filled X: check if it is valid
                        if(check_cross(p,X) == 1)
                            % This ordering is valid, so let's use the
                            % values in the lam vector to actually
                            % construct the matrix of eigenvalues from X.
                            L{ct} = eigen_from_order(X,lam,dim);
                            if((esc_if_npos == 1 && ~IsPSD(L{ct})) || (lim > 0 && ct >= lim))
                                esc_now = 1;
                                return;
                            end
                            ct = ct + 1;
                        end
                    elseif(j == p) % we are at the end of a column: move to the next row
                        fill_matrix(i+1,i+1,3);
                    else % we are somewhere in the middle of the matrix: move to the next column
                        fill_matrix(i,j+1,k+1);
                    end
                    if(esc_now == 1)
                        return;
                    end
                end
                num_pool(k) = 1;
            end
        end
        X(i,j) = p*(p+1)/2+1; % anything larger than p*(p+1)/2
    end
end

% This function makes sure that the entry of X that is currently being
% placed does not violate the fact that the columns of X must be ascending
% (the rows will automatically be ascending by the method of construction
% of X).
function is_o = check_ordered(X,i,j)
    if(i > 1 && X(i-1,j) > X(i,j))
        is_o = 0;
    else
        is_o = 1;
    end
end

% This function makes sure that X does not contain any of the "criss-cross"
% suborderings described in the "A Better Upper Bound" section of
% http://www.njohnston.ca/2014/02/counting-the-possible-orderings-of-pairwise-multiplication/
function cc = check_cross(p,X)
    for i=1:p % i: row
        for j=1:p % j: column
            for k=1:p % k: row
                for l=1:p % l: column
                    for m=1:p % m: row
                        for n=1:p % n: column
                            if(X(min(i,j),max(i,j)) > X(min(k,l),max(k,l)) && X(min(l,n),max(l,n)) > X(min(j,m),max(j,m)) && X(min(i,n),max(i,n)) < X(min(k,m),max(k,m)))
                                cc = 0;
                                return
                            end
                        end
                    end
                end
            end
        end
    end

    cc = 1;
end

% This function converts an ordering matrix to the corresponding eigenvalue
% constraint matrix. For example, the ordering matrix [1 2 3;* 4 5;* * 6]
% corresponds to the following eigenvalue constraint matrix:
% [2*lam(9), lam(8)-lam(1), lam(7)-lam(2);
%    *     , 2*lam(6)     , lam(5)-lam(3);
%    *     ,     *        , 2*lam(4)]
% For details of how these eigenvalue constraint matrices are constructed,
% see the Hildebrand paper.
function L_mat = eigen_from_order(X,lam,dim)
    prod_dim = prod(dim);
    
    % Do some pre-processing to extract information about the ordering X
    % that will be useful for us.
    Y = triu(X,1);
    X = triu(X) + Y.';
    sparsity_pattern = logical(Y);
    
    % Reduce the entries in Y so that they are 1, 2, ..., (p-1)*p/2, but in
    % the same order as the original entries of Y.
    [~,ind] = sort(nonzeros(Y(sparsity_pattern).'));
    Y(sparsity_pattern) = perm_inv(ind); % Y now contains the indices of the eigenvalues that are subtracted in the off-diagonal entries
    Y = Y + Y' + diag(prod_dim + 1 - diag(X));
    
    % We have now built all of the pieces we need: put them together to
    % make the eigenvalue matrix in three steps.
    L_mat = lam(prod_dim + 1 - X); % this places all of the eigenvalues that have a positive coefficient
    L_mat = L_mat + 2*diag(diag(L_mat)); % this triples the diagonal (we only really want to double it, but one copy of it will be subtracted at the next step)
    L_mat = L_mat - lam(Y); % this takes off all of the eigenvalues that have a negative coefficient
end