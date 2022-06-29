%%  SYMIND    Creates all possible non-decreasing integer vectors
%   This function has two required arguments:
%     D: a positive integer (the number of entries in the desired vectors)
%     V: a vector with increasing entries, like 1:N
%
%   SI = symind(D,V) is a matrix whose rows are all possible vectors of the
%   form (V(1),V(2),...,V(D)), where V(1) <= V(2) <= ... <= V(D). These
%   rows are arranged in SI in lexicographical order.
%
%   URL: http://www.qetlab.com/symind

%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: May 10, 2022

function si = symind(d,v)
    n = length(v);
    if(d == 1)
        si = v(:);
    else
        ci = 0;
        si = zeros(nchoosek(n+d-1,d),d);
        for j = 1:n
            ti = symind(d-1,v(v>=v(j)));
            cia = size(ti,1);
            si(ci+1:ci+cia,:) = [v(j)*ones(cia,1),ti];
            ci=ci+cia;
        end
    end
end