%%  COMPOUNDMATRIX    Computes the rth compound matrix of a given matrix
%   This function has two required arguments:
%     A: an arbitrary matrix
%     r: a positive integer
%
%   comp = CompoundMatrix(A, r) returns the rth compound matrix of A
%
%   If r > min{m, n}, then the rth compound matrix of A is the 0x0 matrix.
%   Otherwise, the size of the result is (m choose r) x (n choose r).
%
%   URL: http://www.qetlab.com/CompoundMatrix
%   !! CURRENTLY THERE IS NO PAGE FOR THIS FUNCTION !!
%             
%   author: Benjamin Talbot
%   package: QETLAB
%   last updated: August 8, 2024

function comp = CompoundMatrix(A, r)
    m = size(A, 1);
    n = size(A, 2);
    if r > min(m, n)
        comp = [];
        return
    end

    rows = nchoosek(1:m, r);
    cols = nchoosek(1:n, r);
    dims = [size(rows, 1), size(cols, 1)];
    comp = zeros(dims);
    for i = 1:dims(1)
        for j = 1:dims(2)
            comp(i, j) = det(A(rows(i, :), cols(j, :)));
        end
    end
end
