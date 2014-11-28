%%  ANTISYMMETRICPROJECTION    Produces the projection onto the antisymmetric subspace
%   This function has one required argument:
%     DIM: the dimension of the local systems
%
%   PA = AntisymmetricProjection(DIM) is the orthogonal projection onto the
%   anti symmetric subspace of two copies of DIM-dimensional space. PA is
%   always a sparse matrix.
%
%   This function has five optional arguments:
%     P (default 2)
%     PARTIAL (default 0)
%     MODE (default -1)
%
%   PA = AntisymmetricProjection(DIM,P,PARTIAL,MODE) is the orthogonal
%   projection onto the antisymmetric subspace of P copies of
%   DIM-dimensional space. If PARTIAL = 1 then PA isn't the orthogonal
%   projection itself, but rather a matrix whose columns form an
%   orthonormal basis for the antisymmetric subspace (and hence PA*PA' is
%   the orthogonal projection onto the antisymmetric subspace). MODE is a
%   flag that determines which of two algorithms is used to compute the
%   antisymmetric projection. If MODE = -1 then this script chooses
%   whichever algorithm it thinks will be faster based on the values of DIM
%   and P. If you wish to force one specific algorithm, set either MODE = 0
%   (which generally works best when DIM is small) or MODE = 1 (which
%   generally works best when P is small).
%
%   URL: http://www.qetlab.com/AntisymmetricProjection

%   requires: opt_args.m, perm_sign.m, PermutationOperator.m,
%             PermuteSystems.m, sporth.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: November 26, 2012

function PA = AntisymmetricProjection(dim,varargin)

    % set optional argument defaults: p=2, partial=0, mode=-1
    [p,partial,mode] = opt_args({ 2, 0, -1 },varargin{:});
    dimp = dim^p;
    
    if(p == 1)
        PA = speye(dim);
        return
    elseif(dim < p) % the antisymmetric subspace is empty if dim < p
        PA = sparse(dimp,dimp*(1-partial));
        return
    elseif(mode == -1) % guess which algorithm will be faster
        mode = (dim >= p+4);
    end
    
    plist = perms(1:p);
    pfac = size(plist,1);

    % The antisymmetric projection is the normalization of the sum of all
    % signed permutation operators. This method of computing the
    % antisymmetric projection is reasonably quick when the local dimension
    % is large compared to the number of systems.
    if(mode == 1)
        PA = sparse(dimp,dimp);
        for j = 1:pfac
            PA = PA + perm_sign(plist(j,:))*PermutationOperator(dim*ones(1,p),plist(j,:),0,1);
        end
        PA = PA/pfac;
        
        if(partial)
            PA = sporth(PA);
        end
        
    % When the local dimension is small compared to the number of systems,
    % it is quicker to just explicitly construct an orthonormal basis for
    % the antisymmetric subspace.
    else
        lim = nchoosek(dim,p);
        ind = nchoosek(1:dim,p);
        Id = speye(dim);
        PA = spalloc(dimp,lim,pfac*lim); % allocate the correct number of non-zero elements to PA for speed/memory reasons

        % construct an orthonormal basis for the antisymmetric subspace
        for j = 1:lim
            v = spalloc(dimp,1,pfac);
            for k = 1:pfac
                vt = Id(:,ind(j,plist(k,1)));
                for m = 2:p
                    vt = kron(vt,Id(:,ind(j,plist(k,m))));
                end
                v = v + perm_sign(plist(k,:))*vt;
            end
            PA(:,j) = v;
        end

        if(partial)
            PA = PA/sqrt(pfac);
        else
            PA = PA*PA'/pfac;
        end
    end
end