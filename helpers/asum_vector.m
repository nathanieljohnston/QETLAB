%%  ASUM_VECTOR    Creates all binary vectors with a given sum
%   This function has two required arguments:
%     SM: a positive integer (the desired sum)
%     P: a positive integer (the number of entries in the desired vectors)
%
%   ALIST = asum_vector(SM,P) is a matrix whose rows are all vectors with P
%   0 or 1 entries adding up to SM.
%
%   This function is a companion to asymind, in much the same way that
%   sum_vector is a companion to symind.
%
%   URL: http://www.qetlab.com/asum_vector

%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: June 29, 2022

function alist = asum_vector(sm,p)
    if p <= 1
        if(sm <= 1)
            alist = sm;
        else
            alist = zeros(0,p);
        end
    else
        k = 0;
        alist = zeros(nchoosek(p,sm),p);
        for j = 1:-1:0
            if(p-1 >= sm-j && sm-j >= 0)
                cs = nchoosek(p-1,sm-j);
                t = [j*ones(cs,1),asum_vector(sm-j,p-1)];
                alist(k+1:k+cs,:) = t;
                k = k + cs;
            end
        end
    end
end