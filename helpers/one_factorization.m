%%  ONE_FACTORIZATION    Computes a 1-factorization of a list of objects
%   This function has one required argument:
%     N: either an even natural number (the number of objects to be
%        matched) or a vector containing an even number of distinct objects
%        to be matched
%
%   FAC = one_factorization(N) is a matrix with each row corresponding to a
%   perfect matching of N objects (or, if N is a vector, each row of PM is
%   a pefect matching of the entries of N). Each perfect matching is read
%   "naively": for each j, PM(j,1) is matched with PM(j,2), PM(j,3) is
%   matched with PM(j,4), and so on. The perfect matchings that are
%   returned have the property that every pair of objects is included in
%   exactly one of the perfect matchings (i.e., these perfect matchings
%   together make up a 1-factorization of the complete graph).
%
%   URL: http://www.qetlab.com/one_factorization

%   requires: nothing
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: November 6, 2014

function fac = one_factorization(n)

if(length(n) == 1)
    n = 1:n;
end
sz = length(n);

% If we were given an odd number of objects, no 1-factorization exists.
if(mod(sz,2) == 1)
    error('one_factorization:DoesNotExist','There is no 1-factorization of an odd number of objects.');
end

fac = zeros(sz-1,sz);
for j=1:sz-1
    fac(j,[1,2]) = [n(j),n(sz)];
    for k=2:(sz/2)
        fac(j,[2*k-1,2*k]) = [n(mod(j-(k-1)-1,sz-1)+1),n(mod(j+(k-1)-1,sz-1)+1)];
    end
end