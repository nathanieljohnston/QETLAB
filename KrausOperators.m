%%  KRAUSOPERATORS    Computes a set of Kraus operators for a superoperator
%   This function has one required input argument:
%     PHI: a superoperator
%
%   KO = KrausOperators(PHI) is a cell containing Kraus operators for PHI.
%   This set of Kraus operators will always be canonical in the following
%   ways:
%     (1) If PHI is completely positive, only the left Kraus operators will
%     be returned (and the right Kraus operators are the same).
%     (2) If PHI is Hermiticity preserving, the right Kraus operators are
%     the same as the left Kraus operators, up to sign. The pairs of Kraus
%     operators that are equal are given first, followed by the pairs that
%     are negatives of each other.
%     (3) The left Kraus operators form an orthogonal set in the
%     Hilbert-Schmidt inner product, and similarly for the right Kraus
%     operators.
%
%   This function has one optional input argument:
%     DIM (default has both subsystems of equal dimension)
%
%   KO = KrausOperators(PHI,DIM) is the same as above, where DIM is a
%   1-by-2 vector containing the input and output dimensions of PHI, in
%   that order (equivalently, these are the dimensions of the first and
%   second subsystems of the Choi matrix PHI, in that order). DIM is
%   required if and only if PHI has unequal input and output dimensions and
%   is provided as a Choi matrix.
%
%   See also: ChoiMatrix
%
%   URL: http://www.qetlab.com/KrausOperators

%   requires: ApplyMap.m, ChoiMatrix.m, iden.m, IsCP.m, IsHermPreserving.m,
%             IsPSD.m, MaxEntangled.m, opt_args.m, PermuteSystems.m
%
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   version: 0.50
%   last updated: January 8, 2013

function ko = KrausOperators(Phi,varargin)

if(iscell(Phi)) % is already a set of Kraus operators, but we want a Canonical set
    dim = [size(Phi{1,1},2),size(Phi{1,1},1)];
    if(size(Phi,2) > 1)
        dim = [dim;size(Phi{1,2},2),size(Phi{1,2},1)];
        dX = [dim(1,1)*dim(1,2),dim(2,1)*dim(2,2)];
    else
        dX = prod(dim)*ones(1,2);
    end
    Phi = ChoiMatrix(Phi);
else % Phi was provided as a Choi matrix
    dX = size(Phi);
    sdX = round(sqrt(dX));

    % Set optional argument defaults: dim=sqrt(length(Phi))
    [dim] = opt_args({ [sdX(1) sdX(1);sdX(2) sdX(2)] },varargin{:});
end

% allow the user to enter a single number for dim
if(length(dim) == 1)
    dim = [dim,dX(1)/dim];
    if abs(dim(2) - round(dim(2))) >= 2*dX(1)*eps
        error('KrausOperators:InvalidDim','If DIM is a scalar, PHI must be square and DIM must evenly divide length(PHI); please provide the DIM array containing the dimensions of the subsystems.');
    end
    dim(2) = round(dim(2));
end

% allow the user to enter a vector for dim if X is square
if(min(size(dim)) == 1)
    dim = dim(:).'; % force dim to be a row vector
    dim = [dim;dim];
end

% Compute a canonical set of Kraus operators for Hermiticity-preserving
% maps.
if(IsHermPreserving(Phi))
    ep = max(dX)*eps(norm(Phi,'fro'));
    [V,S] = eig(full(Phi));
    [S,ind] = sort(diag(S),'descend');
    V = V(:,ind);
    
    ind1 = find(S <= ep,1) - 1;
    ind2 = find(S < -ep,1) - 1;
    if(min(size(ind1)) == 0)
        ind1 = length(S);
    end
    if(min(size(ind2)) == 0)
        ind2 = length(S);
    end
    
    sgnS = sign(S);
    V = V*diag(sqrt(abs(S)));
    ko(:,1) = mat2cell(reshape(V(:,[1:ind1,ind2+1:end]),dim(1,2),dim(1,1)*(ind1+length(S)-ind2)),dim(1,2),dim(1,1)*ones((ind1+length(S)-ind2),1)).';
    if(~IsCP(Phi))
        V = V*diag(sgnS);
        ko(:,2) = mat2cell(reshape(V(:,[1:ind1,ind2+1:end]),dim(1,2),dim(1,1)*(ind1+length(S)-ind2)),dim(1,2),dim(1,1)*ones((ind1+length(S)-ind2),1)).';
    end
    
% Compute a canonical set of Kraus operators for all other maps.
else
    [U,S,V] = svd(full(Phi));
    S = diag(S);
    ind = find(S <= max(dX)*eps(norm(Phi,'fro')),1) - 1;
    if(min(size(ind)) == 0)
        ind = length(S);
    end
    S = diag(sqrt(S));
    U = U*padarray(S,size(U,2)-size(S,1),'post');
    V = V*padarray(S,size(V,2)-size(S,1),'post');

    ko(:,1) = mat2cell(reshape(U(:,1:ind),dim(1,2),dim(1,1)*ind),dim(1,2),dim(1,1)*ones(ind,1)).';
    ko(:,2) = mat2cell(reshape(V(:,1:ind),dim(2,2),dim(2,1)*ind),dim(2,2),dim(2,1)*ones(ind,1)).';
end
