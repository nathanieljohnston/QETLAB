%%  VEC_PARTITIONS    Produces all possible partitions of a vector
%
%   PT = vec_partitions(V,P,SZ) computes all partitions of the vector V
%   into P parts. PT is a cell, each of whose entries is a 1-by-P cell
%   containing the P different parts.
%
%   This function has one optional argument:
%     SZ (default [1,1,...,1])
%   
%   PT = vec_partitions(V,P,SZ) computes all partitions of the vector V
%   into P parts, with the restriction that, for all j, the j-th part must
%   have at least SZ(j) elements.
%
%   URL: http://www.qetlab.com/vec_partitions

%   requires: opt_args.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: October 18, 2014

function pt = vec_partitions(v,p,varargin)

% set optional argument defaults: sz = [1,1,...,1]
[sz] = opt_args({ ones(1,p) },varargin{:});

if(p > 1) % if p > 1 then we still have work to do
    pt = cell(0);
    ub = length(v) - sum(sz(2:end));
    for jj = sz(1):ub
        tpt = nchoosek(v,jj);
        for kk = 1:size(tpt,1) % recurse through the other parties
            tpt2 = vec_partitions(setdiff(v,tpt(kk,:)),p-1,sz(2:end));
            for ll = 1:length(tpt2)
                pt{end+1} = {tpt(kk,:), tpt2{ll}{:}};
            end
        end
    end
else % we are on the last party: only one option left
    pt{1}{1} = v;
end