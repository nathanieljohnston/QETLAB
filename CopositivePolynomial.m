%%  COPOSITIVEPOLYNOMIAL    Creates a homogenous polynomial whose non-negativity is equivalent to copositivity of a given matrix
%   This function has one required argument:
%     C: a matrix
%
%   P = CopositivePolynomial(C) is a vector that represents a homogeneous
%       quartic polynomial whose non-negativity is equivalent to
%       copositivity of P. In particular, this polynomial is x.'*C*x, where
%       x = (x1^2, x2^2, ..., xn^2).
%
%   URL: http://www.qetlab.com/CopositivePolynomial

%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: August 1, 2023

function p = CopositivePolynomial(C)
    n = length(C);
    C = (C+C')/2;

    p = zeros(nchoosek(n+3,n-1),1);% pre-load the polynomial coefficient vector with the correct number of entries
    for i = 1:n
        expind = zeros(1,n);
        expind(i) = 4;
        p(exp2ind(expind)) = C(i,i);
        for j = i+1:n
            expind = zeros(1,n);
            expind(i) = 2;
            expind(j) = 2;
            p(exp2ind(expind)) = 2*C(i,j);
        end
    end
end