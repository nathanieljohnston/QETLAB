function test()
m = 15;
n = 15;
% dim = [m, n];
% X = RandomDensityMatrix(prod(dim));
dim = [2, 2];
X = [1, 0, 0, 1; 0, 0, 0, 0; 0, 0, 0, 0; 1, 0, 0, 1];

tic
[s, U, V] = OperatorSchmidtDecomposition(X, dim);
toc
allhermitian = all(cellfun(@isHermitian, U));
if allhermitian
    allhermitian = all(cellfun(@isHermitian, V));
end
allhermitian

tic
[s1, U1, V1] = OperatorSchmidtDecompositionOld(X, dim);
toc
allhermitian = all(cellfun(@isHermitian, U1));
if allhermitian
    allhermitian = all(cellfun(@isHermitian, V1));
end
allhermitian;

res = TensorSum(s, U, V);
res1 = TensorSum(s1, U1, V1);
all(all(abs(res - X) < 1e-10))
all(all(abs(res1 - X) < 1e-10));

function isHermitian = isHermitian(X)
    isHermitian = all(all(abs(X - X') < 1e-10));
end
end