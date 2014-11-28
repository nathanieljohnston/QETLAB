%%  ISPRODUCTVECTOR   Determines if a pure state is a product vector
%   This function has one required argument:
%     VEC: a vector living in the tensor product of two or more subsystems
%
%   IPV = IsProductVector(VEC) is either 1 or 0, indicating that VEC is or
%   is not a product state (note that VEC is assumed to be bipartite unless
%   the optional argument DIM (see below) is specified).
%
%   This function has one optional input argument:
%     DIM (default has two subsystems of equal dimension)
%
%   [IPV,DEC] = IsProductVector(VEC,DIM) indicates that VEC is or is not a
%   product state, as above. DIM is a vector containing the dimensions of
%   the subsystems that VEC lives on. If IPV = 1 then DEC is a product
%   decomposition of VEC. More specifically, DEC is a cell containing
%   vectors whose tensor product equals VEC.
%
%   URL: http://www.qetlab.com/IsProductVector

%   requires: opt_args.m, SchmidtDecomposition.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: November 26, 2012

function [ipv,dec] = IsProductVector(vec,varargin)

lv = length(vec);

% set optional argument defaults: dim=sqrt(length(vec))
[dim] = opt_args({ round(sqrt(lv)) },varargin{:});

% allow the user to enter a single number for dim
num_sys = length(dim);
if(num_sys == 1)
    dim = [dim,lv/dim];
    if abs(dim(2) - round(dim(2))) >= 2*lv*eps
        error('IsProductVector:InvalidDim','The value of DIM must evenly divide length(VEC); please provide a DIM array containing the dimensions of the subsystems.');
    end
    dim(2) = round(dim(2));
    num_sys = 2;
end

% if there are only two subsystems, just use the Schmidt decomposition
if(num_sys == 2)
    [s,u,v] = SchmidtDecomposition(vec,dim,2);
    ipv = (s(2) <= prod(dim) * eps(s(1)));
    if(ipv) % provide this even if not requested, since it is needed if this function was called as part of its recursive algorithm (see below)
        u = u*sqrt(s(1));
        v = v*sqrt(s(1));
        dec = {u(:,1) v(:,1)};
    end
    
% if there are more subsystems, recursively use the Schmidt decomposition across many cuts until we are sure
else
    [ipv,dec] = IsProductVector(vec,[dim(1)*dim(2),dim(3:end)]);
    if(ipv)
        [ipv,tdec] = IsProductVector(dec{1},[dim(1),dim(2)]);
        if(ipv)
            dec = [tdec dec{2:end}];
        end
    end
end

if(~ipv)
    dec = 0;
end