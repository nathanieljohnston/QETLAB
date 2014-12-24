%%  TENSORSUM  Computes a vector or operator from its tensor decomposition
%   This function can be called in two ways:
%   TS = TensorSum(A1,A2,...)
%   TS = TensorSum(S,A1,A2,...)
%
%   Each term Ai is either a matrix or a cell containing matrices (and they
%   should all be the same: either all matrices or all cells). If they are
%   matrices, then the k-th column of each Ai will be tensored together for
%   all k, and then the sum over k will be taken at the end. If they are
%   cells, then the k-th element of each Ai will be tensored together for
%   all k, and then the sum over k will be taken at the end.
%
%   If S is provided, it is a vector of weights (such as Schmidt
%   coefficients) that will be applied when summing the terms at the end of
%   the computation.
%
%   This function acts as an inverse of the SchmidtDecomposition and
%   OperatorSchmidtDecomposition functions. For example, if [S,U,V] =
%   SchmidtDecomposition(VEC) then VEC = TensorSum(S,U,V).
%
%   URL: http://www.qetlab.com/TensorSum

%   requires: nothing
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: December 19, 2014

function ts = TensorSum(a1,a2,varargin)

% Determine the number of terms in the tensor decomposition.
if(iscell(a2))
    len = length(a2);
else
    len = size(a2,2);
end

% Determine whether or not a vector of coefficients was provided, then read
% in the input arguments.
if(iscell(a1) || min(size(a1)) > 1 || len == 1)
    s = ones(len,1);
    A = [{a1} {a2} varargin];
else
    s = a1;
    A = [{a2} varargin];
end

% Get some more preliminary values (dimensions, etc.) and then convert
% everything into cell matrices (instead of matrices of column vectors) if
% necessary.
num_parties = length(A);
dim = zeros(2,num_parties);
for j = 1:num_parties
    if(~iscell(A{j}))
        A{j} = mat2cell(A{j},size(A{j},1),ones(1,size(A{j},2)));
    end
    dim(:,j) = size(A{j}{1}).';
end

% Finally, actually compute the tensor sum.
ts = sparse(prod(dim(1,:)),prod(dim(2,:)));
for j = 1:len
    ts_tmp = A{1}{j};
    for k = 2:num_parties
        ts_tmp = kron(ts_tmp,A{k}{j});
    end
    ts = ts + s(j)*ts_tmp;
end