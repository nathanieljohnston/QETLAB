%%  POLYNOMIALSOS    Bounds the optimal value of a homogeneous polynomial on the unit sphere via the Sum-Of-Squares hierarchy
%   This function has four required arguments:
%     P: the polynomial to optimize, as a vector of its coefficients in
%        lexicographical order
%     N: the number of variables
%     D: half the degree of the polynomial
%     K: a non-negative integer that indicates the level of the hierarchy
%        used to bound the optimal value
%
%   OB = PolynomialSOS(P,N,D,K) is an upper bound on the maximum value of
%        the N-variable degree-D polynomial P. Higher values of K result
%        a better upper bound, at the expense of increased computational
%        cost.
%
%   This function has two optional arguments:
%     OPTTYPE (optional, default 'max'): either 'max' or 'min', indicating
%        whether the polynomial should be maximized or minimized
%     TARGET (optional, default 'none'): either 'none' or a numeric value
%        indicating a value that the computation should exit early if it
%        proves that the optimal value of the polynomial is larger or
%        smaller than.
%
%   [OB,IB] = PolynomialSOS(P,N,D,K,OPTTYPE,TARGET) gives an outer bound
%        (OB) and an inner bound (IB) on the optimum value of the
%        N-variable degree-D polynomial P.
%
%   URL: http://www.qetlab.com/PolynomialSOS

%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: August 4, 2023

function [ob,ib] = PolynomialSOS(p,n,d,k,varargin)

    % set optional argument defaults: OPTTYPE='max', TARGET='none'
    [opttype,target] = opt_args({ 'max', 'none' },varargin{:});
    do_max = strcmpi(opttype,'max');% true means maximize, false means minimize
    has_target = isnumeric(target);% whether or not a target value was specified

    ob_start = tic;
    M = PolynomialAsMatrix(p,n,d,k);
    s = length(M);

    P = SymmetricProjection(n,d+k,1,0);
    cvx_begin sdp quiet
        cvx_precision best;
        variable rho(s,s) hermitian
        if(do_max)
            maximize real(trace(rho*M))
        else
            minimize real(trace(rho*M))
        end
        subject to
            rho >= 0;
            trace(rho) == 1;
            PartialTranspose(P*rho*P',1,n*ones(1,d+k)) == P*rho*P';
    cvx_end
    ob = real(cvx_optval);
    ob_end = toc(ob_start);% time spent performing outer bound calculation

    % If requested, compute inner bounds via error bounds. But only do this
    % if there was not a target set that we have already crossed over
    % (i.e., we are maximizing and found an upper bound below the target or
    % we are minimizing and found a lower bound above the target).
    if(nargout > 1)
        if(has_target && ((do_max && ob <= target) || (~do_max && ob >= target)))
            % Set ib to a useless bound; we hit the target and so do not
            % need this bound, but it was still requested, so we have to
            % return it to avoid errors.
            if(do_max)
                ib = -Inf;
            else
                ib = Inf;
            end
        else
            ib_start = tic;
            ib_end = 0;
            if(k > 0)% error estimate is based on the k = 0 matrix
                M = PolynomialAsMatrix(p,n,d);
            end
    
            if(do_max)
                ib = ob - 4*d*(n-1)*norm(full(M))/(d+k+1);
            else
                ib = ob + 4*d*(n-1)*norm(full(M))/(d+k+1);
            end
    
            si = symind(2*d,1:n);
            ob_end = ob_end / 4;
            while(ib_end < ob_end)% keep computing randomized inner bounds until we have spent at least 25% as much time on this as we did on outer bounds
                new_ib = poly_rand_input(p,n,si);
                if(do_max)
                    ib = max(ib,new_ib);
                else
                    ib = min(ib,new_ib);
                end
                ib_end = toc(ib_start);
            end
        end
    end
end