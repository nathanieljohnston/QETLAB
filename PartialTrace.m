%%  PARTIALTRACE    Computes the partial trace of a matrix
%   This function has one required argument:
%     X: a square matrix
%
%   XPT = PartialTrace(X) is the partial transpose of the matrix X,
%   where it is assumed that length(X) is a perfect squares and both
%   subsystems have equal dimension. The trace is taken over the second
%   subsystem.
%
%   This function has three optional arguments:
%     SYS (default 2)
%     DIM (default has all subsystems of equal dimension)
%     MODE (default -1)
%
%   XPT = PartialTrace(X,SYS,DIM,MODE) gives the partial trace of the matrix X,
%   where the dimensions of the (possibly more than 2) subsystems are given
%   by the vector DIM and the subsystems to take the trace on are given by
%   the scalar or vector SYS. MODE is a flag that determines which of two
%   algorithms is used to compute the partial trace. If MODE = -1 then this
%   script chooses whichever algorithm it thinks will be faster based on
%   the dimensions of the subsystems being traced out and the sparsity of
%   X. If you wish to force one specific algorithm, set either MODE = 0
%   (which generally works best for full or non-numeric matrices, or sparse
%   matrices when most of the subsystems are being traced out) or MODE = 1
%   (which generally works best when X is large and sparse, and the partial
%   trace of X will also be large).
%
%   URL: http://www.qetlab.com/PartialTrace

%   requires: opt_args.m, PermuteSystems.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: November 22, 2012

function Xpt = PartialTrace(X,varargin)

lX = length(X);

% set optional argument defaults: sys=2, dim=round(sqrt(length(X))), mode=-1
[sys,dim,mode] = opt_args({ 2, round(sqrt(lX)), -1 },varargin{:});

num_sys = length(dim);

% allow the user to enter a single number for dim
if(num_sys == 1)
    dim = [dim,lX/dim];
    if abs(dim(2) - round(dim(2))) >= 2*lX*eps
        error('PartialTrace:InvalidDim','If DIM is a scalar, DIM must evenly divide length(X).');
    end
    dim(2) = round(dim(2));
    num_sys = 2;
end

sp = issparse(X);
isnum = isnumeric(X);
prod_dim = prod(dim);
prod_dim_sys = prod(dim(sys));

% Determine which of two computation methods to use (i.e., guess which
% method will be faster).
if(mode == -1)
    mode = (isnum && sp && prod_dim_sys^2 <= prod_dim);
end

% If the matrix is sparse and the amount we are tracing over is smaller
% than the amount remaining, just do the naive thing and manually add up
% the blocks.
if(mode)
    sub_sys_vec = prod_dim*ones(1,prod_dim_sys)/prod_dim_sys;
    perm = [sys,setdiff(1:num_sys,sys)];

    X = mat2cell(PermuteSystems(X,perm,dim), sub_sys_vec, sub_sys_vec); % permute the subsystems so that we just have to do the partial transpose on the first (potentially larger) subsystem
    Xpt = sparse(sub_sys_vec(1),sub_sys_vec(1));
    for j = 1:prod_dim_sys
        Xpt = Xpt + X{j,j};
    end
    
% Otherwise, do a clever trick with mat2cell or reshaping, which is almost always faster.
else
    sub_prod = prod_dim/prod_dim_sys;
    sub_sys_vec = prod_dim*ones(1,sub_prod)/sub_prod;
    perm = [setdiff(1:num_sys,sys),sys];

    Xpt = PermuteSystems(X,perm,dim); % permute the subsystems so that we just have to do the partial trace on the second (potentially larger) subsystem

    if(isnum) % if the input is a numeric matrix, perform the partial trace operation the fastest way we know how
        Xpt = cellfun(@(x) full(trace(x)), mat2cell(Xpt, sub_sys_vec, sub_sys_vec)); % partial trace on second subsystem
        if(sp) % if input was sparse, output should be too
            Xpt = sparse(Xpt);
        end
    else % if the input is not numeric (such as a variable in a semidefinite program), do a slower method that avoids mat2cell (mat2cell doesn't like non-numeric arrays)
        Xpt = reshape(permute(reshape(Xpt,[sub_sys_vec(1),sub_prod,sub_sys_vec(1),sub_prod]),[2,4,1,3]),[sub_prod,sub_prod,sub_sys_vec(1)^2]);
        Xpt = sum(Xpt(:,:,1:sub_sys_vec(1)+1:sub_sys_vec(1)^2),3);
    end
end