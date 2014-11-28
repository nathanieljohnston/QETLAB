%%  SCHMIDTDECOMPOSITION   Computes the Schmidt decomposition of a bipartite vector
%   This function has one required argument:
%     VEC: a bipartite vector to have its Schmidt decomposition computed
%
%   S = SchmidtDecomposition(VEC) is a vector containing the non-zero
%   Schmidt coefficients of the bipartite vector VEC, where the two
%   subsystems are each of size sqrt(length(VEC)).
%
%   This function has two optional input arguments:
%     DIM (default [sqrt(length(VEC)),sqrt(length(VEC))])
%     K (default 0)
%
%   [S,U,V] = SchmidtDecomposition(VEC,DIM,K) gives the Schmidt
%   coefficients S of the vector VEC and the corresponding left and right
%   Schmidt vectors in the matrices U and V. DIM is a 1x2 vector containing
%   the dimensions of the subsystems that VEC lives on. K is a flag that
%   determines how many terms in the Schmidt decomposition should be
%   computed. If K = 0 then all terms with non-zero Schmidt coefficients
%   are computed. If K = -1 then all terms (including zero Schmidt
%   coefficients) are computed. If K > 0 then the K terms with largest
%   Schmidt coefficients are computed.
%
%   If DIM is a scalar instead of a vector, then it is assumed that the
%   first subsystem of size DIM and the second subsystem of size
%   length(VEC)/DIM.
%
%   URL: http://www.qetlab.com/SchmidtDecomposition

%   requires: opt_args.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: December 1, 2012

function [s,u,v] = SchmidtDecomposition(vec,varargin)

lv = length(vec);

% set optional argument defaults: dim=sqrt(length(vec)), k=0
[dim,k] = opt_args({ round(sqrt(lv)), 0 },varargin{:});

% allow the user to enter a single number for dim
if(length(dim) == 1)
    dim = [dim,lv/dim];
    if abs(dim(2) - round(dim(2))) >= 2*lv*eps
        error('SchmidtDecomposition:InvalidDim','The value of DIM must evenly divide length(VEC); please provide a DIM array containing the dimensions of the subsystems.');
    end
    dim(2) = round(dim(2));
end

% Try to guess whether svd or svds will be faster, and then perform the
% appropriate singular value decomposition.
adj = 20 + 1000*(~issparse(vec));

if(k > 0 && k <= ceil(min(dim)/adj)) % just a few Schmidt coefficients
    [v,s,u] = svds(reshape(vec,dim(end:-1:1)),k);
else % lots of Schmidt coefficients
    [v,s,u] = svd(reshape(full(vec),dim(end:-1:1)));
    if(k > 0)
        v = v(:,1:k);
        s = s(:,1:k);
        u = u(:,1:k);
    end
end
s = diag(s);
if(k == 0)
    r = sum(s > max(dim) * eps(s(1)));  % Schmidt rank (use same tolerance as MATLAB's rank function)
    s = s(1:r);   % Schmidt coefficients
    u = u(:,1:r);
    v = v(:,1:r);
end
u = conj(u);