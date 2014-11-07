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
%   version: 0.50
%   last updated: November 28, 2012

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

num_parties = length(A);

% If the objects to be tensored are operators, we don't know any tricks to
% speed it up -- just loop through and tensor them all together as
% requested. We could reshape and then tensor as vectors and then use the
% trick below, but the overhead of reshaping seems to make it not worth it.
if(iscell(A{1}))
    dim = zeros(2,num_parties);
    for j = 1:num_parties
        dim(:,j) = size(A{j}{1}).';
    end
    
    % Finally, compute the tensor sum.
    ts = sparse(prod(dim(1,:)),prod(dim(2,:)));
    for j = 1:len
        ts_tmp = A{1}{j};
        for k = 2:num_parties
            ts_tmp = kron(ts_tmp,A{k}{j});
        end
        ts = ts + s(j)*ts_tmp;
    end

% If the objects to be tensored are vectors, we can be faster: once we get
% to the last two parties to be tensored, we can just perform a matrix
% multiplication trick rather than tensoring and then summing. This makes
% the computation significantly faster if the number of parties is low.
else
    ts_tmp = cell(1,len);
    for j = 1:len % faster to just loop than use bsxfun or other tricks
        ts_tmp{j} = s(j)*A{1}(:,j);
        for k = 2:num_parties-1
            ts_tmp{j} = kron(ts_tmp{j},A{k}(:,j));
        end
    end

    ts = A{num_parties}*cell2mat(ts_tmp).';
    ts = ts(:);
end