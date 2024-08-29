%%  OPERATORSCHMIDTDECOMPOSITION   Computes the operator Schmidt decomposition of a bipartite operator
%   This function has one required argument:
%     X: a bipartite operator
%
%   S = OperatorSchmidtDecomposition(X) is a vector containing the non-zero
%   operator Schmidt coefficients of the bipartite operator X, where the
%   two subsystems are each of size sqrt(length(X)) and X is assumed to be
%   square.
%
%   This function has two optional input arguments:
%     DIM (default has both subsystems of equal size)
%     K (default 0)
%
%   [S,U,V] = OperatorSchmidtDecomposition(X,DIM,K) gives the operator
%   Schmidt coefficients S of the operator X and the corresponding left and
%   right tensor operators in the cells U and V. DIM is a 1x2 vector containing
%   the dimensions of the subsystems that X lives on. If X is non-square,
%   different row and column dimensions can be specified by putting the row
%   dimensions in the first row of DIM and the column dimensions in the 
%   second row of DIM. K is a flag that determines how many terms in the
%   operator Schmidt decomposition should be computed. If K = 0 then all
%   terms with non-zero operator Schmidt coefficients are computed. If
%   K = -1 then all terms (including zero terms) are computed. If K > 0
%   then the K terms with largest operator Schmidt coefficients are
%   computed.
%
%   URL: http://www.qetlab.com/OperatorSchmidtDecomposition

%   requires: opt_args.m, PermuteSystems.m, SchmidtDecomposition.m, Swap.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: November 12, 2014

function [s,U,V] = OperatorSchmidtDecomposition(X,varargin)

dX = size(X);
sdX = round(sqrt(dX));

% set optional argument defaults: dim=sqrt(length(X)), k=0
[dim,k] = opt_args({ [sdX(1) sdX(1);sdX(2) sdX(2)], 0 },varargin{:});

% allow the user to enter a single number for dim
if(length(dim) == 1)
    dim = [dim,dX(1)/dim];
    if abs(dim(2) - round(dim(2))) >= 2*dX(1)*eps
        error('OperatorSchmidtDecomposition:InvalidDim','If DIM is a scalar, X must be square and DIM must evenly divide length(X); please provide the DIM array containing the dimensions of the subsystems.');
    end
    dim(2) = round(dim(2));
end

% allow the user to enter a vector for dim if X is square
if(min(size(dim)) == 1)
    dim = dim(:)'; % force dim to be a row vector
    dim = [dim;dim];
end

% The operator Schmidt decomposition is just the Schmidt decomposition of a
% related vector obtained by moving matrix elements around.
[s,u,v] = SchmidtDecomposition(Swap(reshape(X,prod(prod(dim)),1),[2,3],[dim(2,:),dim(1,:)]),prod(dim),k);

% Now reshape things into the proper output format.
sr = length(s);
U = mat2cell(reshape(u,dim(1,1),dim(2,1)*sr),dim(1,1),dim(2,1)*ones(1,sr));
V = mat2cell(reshape(v,dim(1,2),dim(2,2)*sr),dim(1,2),dim(2,2)*ones(1,sr));

% if the input is hermitian, ensure that the output matrices are hermitian
if isHermitian(X)
    allhermitian = all(cellfun(@isHermitian, U));
    if allhermitian
        allhermitian = all(cellfun(@isHermitian, V));
    end

    if allhermitian
        return
    end

    % Generate hermitian bases for the subsystems
    basis1 = genHermBasis(dim(1));
    basis2 = genHermBasis(dim(2));
    dim1_sq = dim(1)^2;
    dim2_sq = dim(2)^2;

    % represent X as a coordinate vector in the hermitian basis
    [j_indices, i_indices] = ndgrid(1:dim1_sq, 1:dim2_sq);
    kron_basis2 = arrayfun(@(j, i) kron(basis1{j}, basis2{i}), j_indices, i_indices, 'UniformOutput', false);
    
    % Reshape the result into a cell array
    kron_basis2 = reshape(kron_basis2, 1, []);
    coords = X(:)' * cell2mat(cellfun(@(B) B(:), kron_basis2, 'UniformOutput', false));
    reshape(coords,[dim1_sq, dim2_sq]);

    % reshape the coordinates into a matrix, compute svd of X in hermitian basis
    reshaped = reshape(coords,[dim1_sq, dim2_sq]);
    [A, s, B] = svd(reshaped);
    s = diag(s);

    % construct U and V from the SVD, ensuring that they are hermitian
    U = cell(1, length(s));
    for i = 1:length(s)
        U{i} = zeros(dim(1));
        for j = 1:dim1_sq
            U{i} = U{i} + A(j, i) * basis1{j};
        end
    end

    V = cell(1, length(s));
    for i = 1:length(s)
        V{i} = zeros(dim(2));
        for j = 1:dim2_sq
            V{i} = V{i} + B(j, i) * basis2{j};
        end
    end
end

function isHermitian = isHermitian(X)
    isHermitian = all(all(abs(X - X') <= 1e-10));
end

function basis = genHermBasis(dim)
    basis = cell(1, dim^2);
    ct = 1;
    for ct = 1:dim
        elem = sparse(dim, dim);
        elem(ct, ct) = 1;
        basis{ct} = elem;
    end
    ct = dim + 1;

    for row = 1:dim
        for col = (row + 1):dim
            elem = sparse(dim, dim);
            elem(row, col) = 1;
            elem(col, row) = 1;
            elem = elem / sqrt(2);
            basis{ct} = elem;
            ct = ct + 1;
        end
    end

    for row = 1:dim
        for col = (row + 1):dim
            elem = sparse(dim, dim);
            elem(row, col) = 1i;
            elem(col, row) = -1i;
            elem = elem / sqrt(2);
            basis{ct} = elem;
            ct = ct + 1;
        end
    end
end
end