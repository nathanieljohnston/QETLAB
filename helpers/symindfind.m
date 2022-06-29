%%  SYMINDFIND    Finds the location of a vector in symind
%   This function has two required arguments:
%     SV: the vector to search for
%     N: the largest possible value that could appear in an entry of SV
%
%   IND = symindfind(SV,N) is the row index of SV in symind(length(SV),1:N)
%
%   URL: http://www.qetlab.com/symindfind

%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: May 10, 2022

function ind = symindfind(sv,n)
    p = length(sv);
    ind = 1;

    for j = 1:p
        for k = 2:sv(j)
            ind = ind + nchoosek(n+p-j-k,p-j);
        end
    end
end