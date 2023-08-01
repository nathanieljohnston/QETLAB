%%  POLYNOMIALASMATRIX    Creates a compact fully symmetric matrix representation of a polynomial
%   This function has four required arguments:
%     P: the polynomial, represented as a vector of coefficients of its
%        monomials in lexicographic order
%     N: number of variables in the polynomial
%     D: half the degree of the polynomial
%
%   GAMMA = Purity(RHO) is the purity of the quantum state RHO (i.e., GAMMA
%   is the quantity trace(RHO^2)).
%
%   This function has one optional argument:
%     K (default 0): A non-negative integer that indicates the level of the
%                    SOS or SOS-type hierarchy. More specifically, this
%                    input argument causes the output matrix to represent
%                    the polynomial (S(x))^K * P(x) instead of p(x) itself,
%                    where S(x) = x1^2 + x2^2 + ... + xN^2.
%
%   URL: http://www.qetlab.com/PolynomialAsMatrix

%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: July 31, 2023

function M = PolynomialAsMatrix(p,n,d,varargin)

    % set optional argument defaults: K = 0
    [k] = opt_args({ 0 },varargin{:});
    
    k = k+d;
    lk = nchoosek(n+k-1,k);
    rki = zeros(1,n);
    rki(1) = k;

    SId = symind(d,1:n);
    rd = countOcc(SId,n);
    ld = nchoosek(n+d-1,d);

    % pre-compute logs of factorials
    logFac = [0,cumsum(log(1:max(2*d,k)))];% This trick lets us compute these logs of factorials a bit faster than computing each one individually; we re-use a lot of values.

    M = sparse(lk,lk);
    for i = 1:lk
        for j = 1:ld
            rh = rki - rd(j,:);
            if(all(rh >= 0) && all(rh <= k-d))
                for kk = 1:ld
                    pVal = p(symindfind(sort([SId(j,:),SId(kk,:)]),n));
                    if(pVal ~= 0)% Only do these computations if the relevant coefficient of the polynomial is non-zero; otherwise it is a waste of time.
                        rkl = rh + rd(kk,:);
                        l = symindfind(spreadOcc(rkl),n);
                        M(i,l) = M(i,l) + p(symindfind(sort([SId(j,:),SId(kk,:)]),n)) * exp(lnM(rh,k-d,logFac) - lnM(rd(j,:)+rd(kk,:),2*d,logFac) - logFac(k+1) + lnG(rki,rh,d,logFac) + lnG(rkl,rh,d,logFac));
                    end
                end
            end
        end
        rki = nextCVec(rki);
    end
end

% Computes the "r" quantity: vectors describing how often entry j occurs in
% a symmetric index.
function r = countOcc(si,n)
    for j = n:-1:1
        r(:,j) = sum(si==j,2);
    end
end

function nx = nextCVec(ox)
    n = length(ox);
    nx = ox;
    addamt = nx(n) + 1;
    nx(n) = 0;

    for j = n-1:-1:1
        if(nx(j) > 0)
            nx(j) = nx(j) - 1;
            nx(j+1) = nx(j+1) + addamt;
            break;
        end
    end
end

% Inverse of countOcc.
function si = spreadOcc(r)
    k = sum(r);
    n = length(r);
    si = zeros(1,k);

    for j = n:-1:1
        newK = k - r(j);
        si((newK+1):k) = j*ones(1,r(j));
        k = newK;
    end
end

% Computes the "m" function: scaling quantities (multinomial coefficients).
function m = lnM(r,d,logFac)
    m = logFac(d+1) - sum(logFac(r+1));
end

% Computes the "g" function: scaling quantities.
function g = lnG(rI,rH,d,logFac)
    g = logFac(d+1) + sum(logFac(rI+1))/2 - sum(logFac(rI-rH+1));
end