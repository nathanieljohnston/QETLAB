function [Q r] = sporth(S,varargin)
% Q = sporth(S)
% returns an (sparse) orthonormal basis for the range of S.
% That is, Q'*Q = I, the columns of Q span the same space as 
% the columns of S, and the number of columns of Q is the 
% rank of S.
%
% [~, r] = sporth(S)
% returns the rank of S
%
% If S is sparse, Q is obtained from the QR decomposition.
% Otherwise, Q is obtained from the SVD decomposition
%
% Bruno Luong <brunoluong@yahoo.com>
% History
%   10-May-2010: original version
%   13-Dec-2012: added an optional second argument TOL, which specifies the
%   tolerance (Nathaniel Johnston <nathaniel@njohnston.ca>)
%
% See also SPNULL, NULL, QR, SVD, ORTH, RANK

if issparse(S)
    m = size(S,1);
    try
        [Q R E] = qr(S); %#ok %full QR
        if m > 1
            s = diag(R);
        elseif m == 1
            s = R(1);
        else
            s = 0;
        end
        s = abs(s);
        if(nargin == 1)
            tol = norm(S,'fro') * eps(class(S));
        else
            tol = varargin{1};
        end
        r = sum(s > tol);
        Q = Q(:,1:r);
    catch %#ok
        % sparse QR is not available on old Matlab versions
        err = lasterror(); %#ok
        if strcmp(err.identifier, 'MATLAB:maxlhs')
            Q = orth(full(S));
        else
            rethrow(err);
        end
    end
else % Full matrix
    Q = orth(S);
end

r = size(Q,2);

end