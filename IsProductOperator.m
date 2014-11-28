%%  ISPRODUCTOPERATOR   Determines if an operator is an elementary tensor
%   This function has one required argument:
%     X: an operator living on the tensor product of two or more subsystems
%
%   IPO = IsProductOperator(X) is either 1 or 0, indicating that X is or
%   is not a product operator (note that X is assumed to be bipartite
%   unless the optional argument DIM (see below) is specified).
%
%   This function has one optional input argument:
%     DIM (default has two subsystems of equal dimension)
%
%   [IPO,DEC] = IsProductOperator(X,DIM) indicates that X is or is not a
%   product operator, as above. DIM is a vector containing the dimensions
%   of the subsystems that X acts on. If IPO = 1 then DEC is a product
%   decomposition of X. More specifically, DEC is a cell containing
%   operators whose tensor product equals X.
%
%   URL: http://www.qetlab.com/IsProductOperator

%   requires: opt_args.m, IsProductVector.m, PermuteSystems.m,
%             SchmidtDecomposition.m, Swap.m
%
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: November 12, 2014

function [ipo,dec] = IsProductOperator(X,varargin)

dX = size(X);
sdX = round(sqrt(dX));

% set optional argument defaults: dim=sqrt(length(X))
[dim] = opt_args({ [sdX(1) sdX(1);sdX(2) sdX(2)] },varargin{:});

% allow the user to enter a single number for dim
num_sys = length(dim);
if(num_sys == 1)
    dim = [dim,dX(1)/dim];
    if abs(dim(2) - round(dim(2))) >= 2*dX(1)*eps
        error('IsProductOperator:InvalidDim','If DIM is a scalar, X must be square and DIM must evenly divide length(X); please provide the DIM array containing the dimensions of the subsystems.');
    end
    dim(2) = round(dim(2));
    num_sys = 2;
end

% allow the user to enter a vector for dim if X is square
if(min(size(dim)) == 1)
    dim = dim(:)'; % force dim to be a row vector
    dim = [dim;dim];
end

% reshape the operator into the appropriate vector and then test if it's a product vector
[ipo,dec] = IsProductVector(PermuteSystems(reshape(X,prod(prod(dim)),1),Swap(1:2*num_sys,[1,2],[2,num_sys]),[dim(2,:),dim(1,:)]),prod(dim));

% reshape the decomposition into the proper form
if(ipo)
    dec = cellfun(@(x,y) reshape(x,y(1),y(2)),dec,mat2cell(dim,2,ones(1,num_sys)),'un',0);
end