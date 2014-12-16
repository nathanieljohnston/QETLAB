%%  SYMMETRICPROJECTION    Produces the projection onto the symmetric subspace
%   This function has one required argument:
%     DIM: the dimension of the local systems
%
%   PS = SymmetricProjection(DIM) is the orthogonal projection onto the
%   symmetric subspace of two copies of DIM-dimensional space. PS is always
%   a sparse matrix.
%
%   This function has five optional arguments:
%     P (default 2)
%     PARTIAL (default 0)
%     MODE (default -1)
%
%   PS = SymmetricProjection(DIM,P,PARTIAL,MODE) is the orthogonal
%   projection onto the symmetric subspace of P copies of DIM-dimensional
%   space. If PARTIAL = 1 then PS isn't the orthogonal projection itself,
%   but rather a matrix whose columns form an orthonormal basis for the
%   symmetric subspace (and hence PS*PS' is the orthogonal projection onto
%   the symmetric subspace). MODE is a flag that determines which of two
%   algorithms is used to compute the symmetric projection. If MODE = -1
%   then this script chooses whichever algorithm it thinks will be faster
%   based on the values of DIM and P. If you wish to force one specific
%   algorithm, set either MODE = 0 (which generally works best when DIM is
%   small) or MODE = 1 (which generally works best when P is small).
%
%   URL: http://www.qetlab.com/SymmetricProjection

%   requires: opt_args.m, PermutationOperator.m, PermuteSystems.m, sporth.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: December 15, 2014

function PS = SymmetricProjection(dim,varargin)

    % set optional argument defaults: p=2, partial=0, mode=-1
    [p,partial,mode] = opt_args({ 2, 0, -1 },varargin{:});
    dimp = dim^p;
    
    if(p == 1)
        PS = speye(dim);
        return
    elseif(mode == -1)
        mode = (dim >= p-1);
    end

    % The symmetric projection is the normalization of the sum of all
    % permutation operators. This method of computing the symmetric
    % projection is reasonably quick when the local dimension is large
    % compared to the number of spaces.
    if(mode == 1)
        plist = perms(1:p);
        pfac = factorial(p);
        PS = sparse(dimp,dimp);

        for j = 1:pfac
            PS = PS + PermutationOperator(dim*ones(1,p),plist(j,:),0,1);
        end
        PS = PS/pfac;
        
        if(partial)
            PS = sporth(PS);
        end
        
    % When the local dimension is small compared to the number of spaces,
    % it is quicker to just explicitly construct an orthonormal basis for
    % the symmetric subspace.
    else
        lim = nchoosek(dim+p-1,dim-1);
        PS = spalloc(dimp,lim,dimp); % allocate dim^p non-zero elements to PS for speed/memory reasons
        slist = sum_vector(p,dim);
        Id = speye(dim);

        for j = 1:lim
            ind = cell2mat(arrayfun(@(x,y) repmat(y,1,x), slist(j,:),1:dim,'un',0));
            plist = unique_perms(ind);
            v = spalloc(dimp,1,nchoosek(p,floor(p/2)));
            sp = size(plist,1);
            for k = 1:sp
                vt = Id(:,plist(k,1));
                for m = 2:p
                    vt = kron(vt,Id(:,plist(k,m)));
                end
                v = v + vt/sqrt(sp);
            end
            PS(:,j) = v;
        end
            
        if(~partial)
            PS = PS*PS';
        end
    end
end

% We need some helper functions to help us through the MODE = 0 algorithm.
function slist = sum_vector(dim,p)
    if p <= 1
        slist = dim;
    else
        k = 0;
        slist = zeros(nchoosek(dim+p-1,p-1),p);
        for j = 0:dim
            cs = nchoosek(j+p-2,p-2);
            t = [sum_vector(j,p-1),(dim-j)*ones(cs,1)];
            slist(k+1:k+cs,:) = t;
            k = k + cs;
        end
    end
end