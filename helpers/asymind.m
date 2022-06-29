%%  ASYMIND    Creates all possible increasing integer vectors
%   This function has two required arguments:
%     D: a positive integer (the number of entries in the desired vectors)
%     V: a vector with increasing entries, like 1:N
%
%   AI = asymind(D,V) is a matrix whose rows are all possible vectors of
%   the form (V(1),V(2),...,V(D)), where V(1) < V(2) < ... < V(D). These
%   rows are arranged in AI in lexicographical order.
%
%   URL: http://www.qetlab.com/asymind

%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: May 10, 2022

function ai = asymind(d,v)
    n = length(v);
    if(d == 1)
        ai = v(:);
    elseif(n >= d)
        ci = 0;
        ai = zeros(nchoosek(n,d),d);
        for j = 1:n-d+1
            ti = asymind(d-1,v(v>v(j)));
            cia = size(ti,1);
            ai(ci+1:ci+cia,:) = [v(j)*ones(cia,1),ti];
            ci=ci+cia;
        end
    else
        ai = zeros(0,n);
    end
end