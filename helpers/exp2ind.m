%%  EXP2IND    Looks up a monomial's lexicographical index based on a list of exponents
%   This function has one required argument:
%     EXPNTS: a vector containing exponents of a monomial in a homogeneous
%             polynomial
%
%   IND = exp2ind(EXPNTS) is the lexicographical index of the term
%         x1^EXPNTS(1) * ... * xn^EXPNTS(n). This is typically used to help
%         create a vector representation of a homogeneous polynomial.
%
%   URL: http://www.qetlab.com/exp2ind

%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: August 1, 2023

function ind = exp2ind(expnts)
    n = length(expnts);
    s = sum(expnts);
    sv = zeros(1,s);

    ct = 1;
    for j = 1:n
        sv(ct:ct+expnts(j)-1) = j*ones(1,expnts(j));
        ct = ct + expnts(j);
    end
    ind = symindfind(sv,n);
end