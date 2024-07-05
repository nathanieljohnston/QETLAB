% if r > min{m, n}, then Cr(A) is the unique 0x0 matrix
% otherwise, has size (m; r) x (n; r) (m choose r x n choose r)?

function compMat = compoundMatrix(A, r)
    m = size(A, 1);
    n = size(A, 2);
    if r > min(m, n)
        compMat = [];
        return
    end

    rows = nchoosek(1:m, r);
    cols = nchoosek(1:n, r);
    dims = [size(rows, 1), size(cols, 1)];
    compMat = zeros(dims);
    for i = 1:dims(1)
        for j = 1:dims(2)
            compMat(i, j) = det(A(rows(i, :), cols(j, :)));
        end
    end
end
