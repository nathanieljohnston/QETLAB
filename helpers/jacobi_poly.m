%%  JACOBI_POLY  Computes the coefficients of Jacobi polynomials
%   This function has three required arguments:
%     A: a real parameter (sometimes called alpha) of the Jacobi polynomials
%     B: a real parameter (sometimes called beta) of the Jacobi polynomials
%     N: the degree of the Jacobi polynomial (a non-negative integer)
%
%   JP = jacobi_poly(A,B,N) is an (N+1)-vector containing the coefficients
%   (starting with the coefficient of the highest-order term) of the Jacobi
%   polynomial P_N^{A,B}.
%
%   URL: http://www.qetlab.com/jacobi_poly

%   requires: nothing
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: March 6, 2013

function jp = jacobi_poly(a,b,n)

    if(n <= 0)
        jp = 1;
    elseif(n <= 1)
        jp = [a+b+2, a-b]/2;
    else
        prev_jp = jacobi_poly(a,b,n-1);
        jp = ((2*n+a+b-1)*(a^2-b^2)*[0,prev_jp] + (2*n+a+b-1)*(2*n+a+b)*(2*n+a+b-2)*[prev_jp,0] - 2*(n+a-1)*(n+b-1)*(2*n+a+b)*[0,0,jacobi_poly(a,b,n-2)]) / (2*n*(n+a+b)*(2*n+a+b-2));
    end
end