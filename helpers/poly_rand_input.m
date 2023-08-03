%%  POLY_RAND_INPUT    Evaluates a homogeneous polynomial on a random input from the unit sphere
%   This function has one required argument:
%     P: a polynomial, as a vector of its coefficients in lexicographical
%        order
%     N: the number of variables
%     SI: the matrix computed by symind(D,1:N), where D is the degree of
%         the polynomial. This should be provided as an input, rather
%         than computed by this function, since this function will
%         typically be called numerous times and so SI should be
%         pre-computed to save time.
%
%   PX = poly_rand_input(P,N,SI) is the value of P(X), where X is a
%        randomly generated point from the unit sphere (according to
%        uniform spherical measure).
%
%   URL: http://www.qetlab.com/poly_rand_input

%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: August 2, 2023

function px = poly_rand_input(p,n,si)
    x = randn(n,1);
    x = x/norm(x);
    xd = prod(x(si),2); % compute all distinct D-fold products of entries of x, in lexicographical order of their subscripts

    px = p.'*xd;
end