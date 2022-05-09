%%  SUM_VECTOR    Creates all non-negative integer vectors with a given sum
%   This function has two required arguments:
%     SM: a positive integer (the desired sum)
%     P: a positive integer (the number of entries in the desired vectors)
%
%   SLIST = sum_vector(SM,P) is a matrix whose rows are all vectors with P
%   non-negative integer entries adding up to SM.
%
%   URL: http://www.qetlab.com/sum_vector

%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: May 9, 2022

function slist = sum_vector(sm,p)
    if p <= 1
        slist = sm;
    else
        k = 0;
        slist = zeros(nchoosek(sm+p-1,p-1),p);
        for j = 0:sm
            cs = nchoosek(j+p-2,p-2);
            t = [(sm-j)*ones(cs,1),sum_vector(j,p-1)];
            slist(k+1:k+cs,:) = t;
            k = k + cs;
        end
    end
end