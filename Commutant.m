%%  COMMUTANT    Computes the commutant of a set of matrices
%   This function has one required argument:
%     A: a matrix, or a cell containing one or more matrices of the same size
%
%   C = Commutant(A) is a cell containing an orthonormal basis (in the
%   Hilbert-Schmidt inner product) for the algebra of matrices that commute
%   with each matrix in the cell A. The elements of C are sparse if and
%   only if the elements of A are sparse.
%
%   URL: http://www.qetlab.com/Commutant

%   requires: opt_args.m, PartialTranspose.m, PermuteSystems.m, spnull.m,
%             Swap.m
%
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: November 12, 2014

function C = Commutant(A)

if(iscell(A)) % allow A to just be a matrix instead of a cell if they only want the commutant of one matrix
    A = [A{:}].';
else
    A = A.';
end
dim = size(A,2);
num_ops = length(A)/dim;

% a sneaky (and fast) way of constructing the commutant that works by
% noting that AX = XA if and only if (kron(A,I) - kron(I,A^T))x = 0, where
% x is the vectorization of X
C = reshape(spnull( kron(PartialTranspose(A,2,[num_ops,dim;1,dim]),speye(dim)) - Swap(kron(speye(dim),A),[1,2],[dim,num_ops,dim],1) ),dim,[]);

% from here on, we're just reshaping the resulting data into matrices of
% the proper size and sparsity
num_comm_ops = size(C,2)/dim;
if(~issparse(A))
    C = full(C);
end
C = mat2cell(PartialTranspose(C,2,[1,dim;num_comm_ops,dim]),dim,dim*ones(1,num_comm_ops));