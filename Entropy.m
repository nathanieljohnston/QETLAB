%%  ENTROPY    Computes the von Neumann entropy of a density matrix
%   This function has one required argument:
%     RHO: a density matrix
%
%   ENT = Entropy(RHO) is the (base 2) von Neumann entropy of RHO.
%
%   This function has one optional input argument:
%     BASE (default 2)
%
%   ENT = Entropy(RHO,BASE) is the von Neumann entropy of RHO, computed
%   with logarithms in the base specified by BASE.
%
%   URL: http://www.qetlab.com/Entropy

%   requires: nothing
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   version: 0.50
%   last updated: September 9, 2014

function ent = Entropy(rho,varargin)

% set optional argument defaults: base=2
[base] = opt_args({ 2 },varargin{:});

lam = eig(rho);
if(base == 2)
    ent = -sum(real(lam.*log2(lam)));
else
    ent = -sum(real(lam.*log(lam)))/log(base);
end