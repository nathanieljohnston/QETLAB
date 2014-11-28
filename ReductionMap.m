%%  REDUCTIONMAP    Produces the reduction map
%   This function has one required argument:
%     DIM: a positive integer (the dimension of the reduction map)
%
%   R = ReductionMap(DIM) is the Choi matrix of the reduction map, which is
%   a positive map on DIM-by-DIM matrices.
%
%   This function has one optional argument:
%     K (default 1)
%   
%   R = ReductionMap(DIM,K) is the Choi matrix of the map defined by
%   R(X) = K*trace(X)*eye(DIM^2) - X. This map is K-positive.
%
%   URL: http://www.qetlab.com/ReductionMap

%   requires: iden.m, MaxEntangled.m, opt_args.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: September 29, 2014

function R = ReductionMap(dim,varargin)

% set optional argument defaults: k=1 (the usual reduction map)
[k] = opt_args({ 1 },varargin{:});

psi = MaxEntangled(dim,1,0);
R = k*speye(dim^2) - psi*psi';