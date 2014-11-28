%%  GENGELLMANN    Produces a generalized Gell-Mann operator
%   This function has three required arguments:
%     IND1 (a nonnegative integer from 0 to DIM-1 inclusive)
%     IND2 (a nonnegative integer from 0 to DIM-1 inclusive)
%     DIM (a positive integer indicating the dimension)
%
%   G = GenGellMann(IND1,IND2,DIM) is a DIM-by-DIM Hermitian operator.
%   These matrices span the entire space of DIM-by-DIM matrices as IND1 and
%   IND2 range from 0 to DIM-1, inclusive, and they generalize the Pauli
%   operators when DIM = 2 and the Gell-Mann operators when DIM = 3.
%
%   This function has one optional argument:
%     SP (default 0)
%
%   G = GenGellMann(IND1,IND2,DIM,SP) is as above, with sparsity of the
%   output determined by the value of SP. If SP = 0 then the output will be
%   full, if SP = 1 then the output will be sparse.
%
%   URL: http://www.qetlab.com/GenGellMann

%   requires: opt_args.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: December 18, 2013

function g = GenGellMann(ind1,ind2,dim,varargin)

% set optional argument defaults: sp=0
[sp] = opt_args({ 0 },varargin{:});

if(ind1 == ind2)
    if(ind1 == 0)
        g = speye(dim);
    else
        g = sqrt(2/(ind1*(ind1+1)))*spdiags([ones(ind1,1);-ind1;zeros(dim-ind1-1,1)],0,dim,dim);
    end
else
    E = sparse(dim,dim);
    E(ind1+1,ind2+1) = 1;
    if(ind1 < ind2)
        g = E + E';
    else
        g = 1i*E - 1i*E';
    end
end

if(~sp)
    g = full(g);
end