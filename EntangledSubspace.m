%%  ENTANGLEDSUBSPACE    Produces a basis of an r-entangled subspace
%   This function has two required arguments:
%     DIM (positive integer indicating the subspace dimension)
%     LOCALDIM (a positive integer or 1-by-2 vector indicating the local dimension(s))
%
%   E = EntangledSubspace(DIM,LOCALDIM) is a matrix whose columns form a
%   basis of an entangled DIM-dimensional subspace of
%   C^(LOCALDIM(1)) \otimes C^(LOCALDIM(2)). If LOCALDIM is a scalar, then
%   it is assumed that both local dimensions are equal to LOCALDIM.
%
%   This function has one optional argument:
%     R (default 1)
%
%   E = EntangledSubspace(DIM,LOCALDIM,R) is as above, but the subspace is
%   R-entangled instead of just entangled (i.e., every member of the
%   subspace has Schmidt rank R+1 or larger).
%
%   URL: http://www.qetlab.com/EntangledSubspace

%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: May 6, 2022

function E = EntangledSubspace(dim,localDim,varargin)

% set optional argument defaults: r=1
[r] = opt_args({ 1 },varargin{:});

% make localDim a vector if it was provided as a scalar
if(length(localDim) == 1)
    localDim = [localDim,localDim];
end

% entangled subspaces only exist in some dimensions
if(dim > (localDim(1) - r)*(localDim(2) - r))
    error('EntangledSubspace:InvalidSize',['No r-entangled subspace of this dimension exists.']);
end

m = min(localDim);
plDim = prod(localDim);
V = fliplr(vander(1:m)); % entries of basis vectors come from Vandermonde matrices
E = sparse(plDim,dim);

% Now construct E. It would be slightly (but only slightly) faster to swap
% the order of these two loops, but we do it this way since it leads to
% slightly better numerical stability and smallest entries in E when dim is
% strictly less than (localDim(1) - r)*(localDim(2) - r).
ct = 1;
for k = 1:m-r % loop through the columns of the Vandermonde matrix that we are putting on this diagonal
    for j = (r+1-localDim(2)):(localDim(1)-r-1) % loop through the relevant diagonals of a matrix, to put Vandermonde columns on them
        ell = length(diag(sparse(localDim(2),localDim(1)),j)); % number of entries in this particular diagonal
        if(k <= ell-r)
            D = V(1:ell,k);
            T = sparse(localDim(2),localDim(1)); % need to be really careful with sizes when local dimensions are unequal (could do these 5 lines in just 1 line if dimensions are equal)
            TD = sparse(diag(D,j)); % put those entries on the correct diagonal
            md1 = min(size(TD,1),localDim(2));
            md2 = min(size(TD,2),localDim(1));
            T(1:md1,1:md2) = T(1:md1,1:md2) + TD(1:md1,1:md2);

            E(:,ct) = reshape(T,plDim,1);
            
            ct = ct + 1;
            if(ct > dim) % we have now found enough columns -- we are done
                return;
            end
        end
    end
end