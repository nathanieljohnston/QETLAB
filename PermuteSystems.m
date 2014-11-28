%%  PERMUTESYSTEMS    Permutes subsystems within a state or operator
%   This function has two required arguments:
%     X: the vector or matrix
%     PERM: a permutation vector
%
%   PX = PermuteSystems(X,PERM) permutes the order of the subsystems
%   of the vector or matrix X according to the permutation vector PERM,
%   where each subsystem is assumed to have equal dimension
%
%   This function has three optional arguments:
%     DIM (default has all subsystems of equal dimension)
%     ROW_ONLY (default 0)
%     INV_PERM (default 0)
%   
%   PX = PermuteSystems(X,PERM,DIM,ROW_ONLY,INV_PERM) permutes the order of
%   the subsystems of the vector or matrix X according to the permutation
%   vector PERM, where the dimensions of the subsystems are given by the
%   vector DIM. If X is non-square and not a vector, different row and
%   column dimensions can be specified by putting the row dimensions in the
%   first row of DIM and the column dimensions in the second row of DIM.
%   If ROW_ONLY is set to 1, then only the rows of X are permuted, but not
%   the columns -- this is equivalent tomultiplying X on the left by the
%   corresponding permutation operator, but not on the right. If ROW_ONLY
%   is set to 1 then DIM only needs to contain the row dimensions of the
%   subsystems, even if X is not square. If INV_PERM equals 1 then the
%   inverse permutation of PERM is applied instead of PERM itself.
%
%   URL: http://www.qetlab.com/PermuteSystems

%   requires: opt_args.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: November 6, 2014

function PX = PermuteSystems(X,perm,varargin)

dX = size(X);
is_vec = (min(dX) == 1);
num_sys = length(perm);
if(is_vec)
    vec_orien = 3 - find(dX == 1, 1); % 1 if column vector, 2 if row vector
end

% set optional argument defaults: dim=round(lX^(1/num_sys)), row_only=0, inv_perm=0
[dim,row_only,inv_perm] = opt_args({ [round(dX(1)^(1/num_sys))*ones(1,num_sys); round(dX(2)^(1/num_sys))*ones(1,num_sys)], 0, 0 },varargin{:});

% allow the user to enter a vector for dim if X is square
if(min(size(dim)) == 1)
    dim_tmp = dim(:)'; % force dim to be a row vector
    if(is_vec)
        dim = ones(2,length(dim));
        dim(vec_orien,:) = dim_tmp;
    else
        dim = [dim_tmp;dim_tmp];
    end
end

prod_dimR = prod(dim(1,:));
prod_dimC = prod(dim(2,:));

% Do some basic input checking.
if length(perm) ~= num_sys
    error('PermuteSystems:InvalidPerm','length(PERM) must equal length(DIM).')
elseif ~all(sort(perm) == 1:num_sys)
    error('PermuteSystems:InvalidPerm','PERM must be a permutation vector.')
elseif(dX(1) ~= prod_dimR || (~row_only && dX(2) ~= prod_dimC))
    error('PermuteSystems:InvalidDim','The dimensions specified in DIM do not agree with the size of X.')
end

% Permuting systems for pure states is easy enough, so just make the vector
% full and then perform the permutation (new-ish versions of MATLAB don't
% like sparse multidimensional arrays).
if(is_vec)
    if(inv_perm)
        PX = reshape(ipermute(reshape(full(X),dim(vec_orien,end:-1:1)),num_sys+1-perm(end:-1:1)),dX);
    else
        PX = reshape(permute(reshape(full(X),dim(vec_orien,end:-1:1)),num_sys+1-perm(end:-1:1)),dX);
    end
    
    % Preserve the sparsity of X.
    if(issparse(X))
        PX = sparse(PX);
    end
    return
end

% If X is not a pure state, it's slightly trickier... do *not* just use the
% same pure state trick with repeated indices though, since that has an
% intermediate step of making the matrix a multidimensional array, which
% you can't do with sparse matrices in new-ish version of MATLAB. The trick
% used here reduces the problem to the pure state version of the problem in
% another way that plays nicely with both full and sparse matrices
row_perm = PermuteSystems(1:dX(1),perm,dim(1,:),0,inv_perm);
PX = X(row_perm,:);
if ~row_only
    col_perm = PermuteSystems(1:dX(2),perm,dim(2,:),0,inv_perm);
    PX = PX(:,col_perm);
end