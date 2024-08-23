%%  ELEMSYMPOLY    Computes the kth elementary symmetric polynomial of a vector.
%   This function has two required arguments:
%     x: any vector of numbers
%     k: a nonnegative integer less than or equal to n, the length of x, which 
%        specifies the degree of the terms in the elementary symmetric polynomial
%
%   res = ElemSymPoly(x) returns the computed value of the 
%   kth elementary symmetric polynomial of the vector x.
%
%   URL: http://www.qetlab.com/ElemSymPoly
%             
%   author: Benjamin Talbot
%   package: QETLAB
%   last updated: August 23, 2024

function res = ElemSymPoly(x, k)
    n = length(x);
    if k > n || k < 0
        error('ElemSymPoly:InvalidInput','k must be a nonnegative integer less than or equal to the length of x.');
    end
    res = 0;
    if k == 0
        res = 1;
    elseif k == 1
        res = sum(x);
    elseif k == n
        res = prod(x);
    else
        % for i = 1:n
        %     res = res + x(i) * ElemSymPoly(x(i+1:n), k-1);
        % end
        idxs = nchoosek(1:n, k);
        for i = 1:size(idxs, 1)
            res = res + prod(x(idxs(i, :)));
        end
    end
end