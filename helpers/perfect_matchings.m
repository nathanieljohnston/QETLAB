%%  PERFECT_MATCHINGS    Gives all perfect matchings of N objects
%   This function has one required argument:
%     N: either an even natural number (the number of objects to be
%        matched) or a vector containing an even number of distinct objects
%        to be matched
%
%   PM = perfect_matchings(N) is a matrix with each row corresponding to a
%   perfect matching of N objects (or, if N is a vector, each row of PM is
%   a pefect matching of the entries of N). Each perfect matching is read
%   "naively": for each j, PM(j,1) is matched with PM(j,2), PM(j,3) is
%   matched with PM(j,4), and so on.
%
%   URL: http://www.qetlab.com/perfect_matchings

%   requires: nothing
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: November 6, 2014

function pm = perfect_matchings(n)

if(length(n) == 1)
    n = 1:n;
end
sz = length(n);

% Base case, n = 2: only one perfect matching.
if(sz == 2)
    pm = n;
    return;
    
% There are no perfect matchings of an odd number of objects.
elseif(mod(sz,2) == 1)
    pm = zeros(0,sz);
    return;
end

% Recursive step: build perfect matchings from smaller ones.

% Only do the recursive step once instead of n-1 times: we will then tweak
% the output n-1 times.
lower_fac = perfect_matchings(n(3:end));
lfac_size = size(lower_fac,1);
pm = zeros(0,sz);

% Now build the perfect matchings we actually want.
for j = 2:sz
    tlower_fac = lower_fac;
    tlower_fac(tlower_fac==n(j)) = n(2);
    pm = [pm;[ones(lfac_size,1)*[n(1),n(j)],tlower_fac]];
end