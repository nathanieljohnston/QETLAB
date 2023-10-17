%%  POLYNOMIALOPTIMIZE    Bounds the optimal value of a homogeneous polynomial on the unit sphere
%   This function has four required arguments:
%     P: the polynomial to optimize, as a vector of its coefficients in
%        lexicographical order
%     N: the number of variables
%     D: half the degree of the polynomial
%     K: a non-negative integer that indicates the level of the hierarchy
%        used to bound the optimal value
%
%   OB = PolynomialOptimize(P,N,D,K) is an upper bound on the maximum value
%        of the N-variable degree-D polynomial P. Higher values of K
%        result a better upper bound, at the expense of increased
%        computational cost.
%
%   This function has two optional arguments:
%     OPTTYPE (optional, default 'max'): either 'max' or 'min', indicating
%        whether the polynomial should be maximized or minimized
%     TARGET (optional, default 'none'): either 'none' or a numeric value
%        indicating a value that the computation should exit early if it
%        proves that the optimal value of the polynomial is larger or
%        smaller than.
%
%   [OB,IB] = PolynomialOptimize(P,N,D,K,OPTTYPE,TARGET) gives an outer
%        bound (OB) and an inner bound (IB) on the optimum value of the
%        N-variable degree-D polynomial P.
%
%   URL: http://www.qetlab.com/PolynomialOptimize

%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: September 2, 2023

function [ob,ib] = PolynomialOptimize(p,n,d,k,varargin)

    % set optional argument defaults: OPTTYPE='max', TARGET='none'
    [opttype,target] = opt_args({ 'max', 'none' },varargin{:});
    has_target = isnumeric(target);% whether or not a target value was specified

    % Just do the minimization version of this optimization.
    % If a maximization was requested, convert it and re-call the function.
    if(strcmpi(opttype,'max'))
        if(nargout > 1)
            [ob,ib] = PolynomialOptimize(-p,n,d,k,'min',target);
            ob = -ob;
            ib = -ib;
        else
            ob = -PolynomialOptimize(-p,n,d,k,'min',target);
        end
        return
    end

    if(nargout > 1)
        ob_start = tic;
    end

    M = PolynomialAsMatrix(p,n,d,k);

    % Compute the polynomial (as a vector of coefficients) (s(x))^d, since
    % we will need that.
    slist = sum_vector(2*d,n);
    plen = length(p);
    ps = zeros(plen,1);% initialize ps (the polynomial (s(x))^d) to have the correct size
    logFac = [0,cumsum(log(1:d))];% logarithms of the factorials from 1 to d
    for j = 1:plen
	    if(all(mod(slist(j,:),2)==0))
		    ps(exp2ind(slist(j,:))) = exp(logFac(d+1) - sum(logFac(1+round(slist(j,:)/2))));
	    end
    end
    MS = PolynomialAsMatrix(ps,n,d,k);

    % Pre-computations are done; get the minimum eigenvalue now!
    if(isa(p,'cvx'))
        MSsqrt = sqrtm(inv(full(MS)));
        ob = lambda_min(MSsqrt*M*MSsqrt);% slower CVX-safe way of getting minimum eigenvalue
        return
    end

    % If the input is not a CVX variable, we can compute this minimum
    % eigenvalue much quicker as follows.
    s = length(M);
    sub_dim = 20 + floor(k/10);% Krylov subspace dimension used in sparse eigenvalue computation. If eigenvalues fail to converge, try increasing this value.
    if(s <= max(sub_dim,500))
        ob = min(real(eig(full(M),full(MS))));
    else
        ob = real(eigs(M,MS,1,'smallestreal','SubspaceDimension',sub_dim));
    end

    % If requested, compute inner bounds via randomization. But only do
    % this if there was not a target set that we have already crossed over
    % (i.e., we are maximizing and found an upper bound below the target or
    % we are minimizing and found a lower bound above the target).
    if(nargout > 1)
        if(has_target && ob >= target)
            % Set ib to a useless bound; we hit the target and so do not
            % need this bound, but it was still requested, so we have to
            % return it to avoid errors.
            ib = Inf;
        else
            ob_end = toc(ob_start);% time spent performing outer bound calculation

            ib_start = tic;
            ib_end = -1;% just to make sure the upcoming while loop starts
            ib = Inf;
    
            si = symind(2*d,1:n);
            ob_end = ob_end / 4;
            while(ib_end < ob_end)% keep computing randomized inner bounds until we have spent at least 25% as much time on this as we did on outer bounds
                new_ib = poly_rand_input(p,n,si);
                ib = min(ib,new_ib);
                ib_end = toc(ib_start);
            end
        end
    end
end