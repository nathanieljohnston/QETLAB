%%  TENSOR   Kronecker tensor product of two or more matrices or vectors
%   This function can be called in one of three ways:
%
%   KRN = Tensor(A,B,C,...) returns the Kronecker product A \otimes B
%   \otimes C \otimes ...
%
%   KRN = Tensor(A,m) returns the Kronecker product of A with itself m
%   times.
%
%   KRN = Tensor(A), if A is a cell, returns the Kronecker product of
%   all matrices contained within A.
%
%   URL: http://www.qetlab.com/Tensor

%   requires: nothing
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: November 28, 2012

function krn = Tensor(A,varargin)

% If A is a cell, tensor together each element of the cell.
if(iscell(A))
    krn = A{1};
    for j = 2:length(A);
        krn = kron(krn,A{j});
    end
    
% If two arguments were received and the second is a scalar, tensor A with itself that many times.
elseif(nargin == 2 && length(varargin{1}) == 1)
    % Tensor naively if we only want a few copies.
    if(varargin{1} <= 4 || (length(A) > 3 && ~issparse(A)))
        krn = A;
        for j = 2:varargin{1};
            krn = kron(krn,A);
        end
        
    % Be more clever (a-la exponentiation by squaring) if we want lots of
    % copies. For n copies, this procedure only does O(log(n)) Kronecker
    % products. Doesn't help much for full matrices, but in practice is
    % about a 50% speedup for some sparse matrices.
    else
        l2 = floor(log2(varargin{1}));
        krn_cell = cell(1,l2+1);
        b = de2bi(varargin{1});
        
        krn = 1;
        krn_cell{1} = A;
        for j = 1:l2
            krn_cell{j+1} = kron(krn_cell{j},krn_cell{j});
            if(b(j) == 1)
                krn = kron(krn,krn_cell{j});
            end
        end
        krn = kron(krn,krn_cell{l2+1});
    end

% If two or more arguments were received and the second isn't a scalar, tensor together all arguments.
else
    krn = A;
    for j = 2:nargin;
        krn = kron(krn,varargin{j-1});
    end
end