%%  APPLYMAP    Applies a superoperator to an operator
%   This function has two required arguments:
%     X: a matrix
%     PHI: a superoperator
%
%   PHIX = ApplyMap(X,PHI) is the operator PHI(X). That is, it is the
%   result of applying the superoperator PHI to the operator X. PHI should
%   be provided either as a Choi matrix, or as a cell with either 1 or 2
%   columns whose entries are its Kraus operators (see full QETLAB
%   documentation for details).
%
%   This function has no optional arguments.
%
%   URL: http://www.qetlab.com/ApplyMap

%   requires: opt_args.m, PermuteSystems.m, Swap.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: November 12, 2014

function PhiX = ApplyMap(X,Phi)

% Both of the following methods of applying the superoperator are much
% faster than naively looping through Kraus operators or constructing
% eigenvectors of a Choi matrix.
if(iscell(Phi)) % the superoperator was given as a cell of Kraus operators
    sPhi = size(Phi);
    if(sPhi(2) == 1 || (sPhi(1) == 1 && sPhi(2) > 2)) % map is CP
        Phi = Phi(:);
        Phi(:,2) = cellfun(@ctranspose,Phi(:,1),'UniformOutput',false);
    else
        Phi(:,2) = cellfun(@ctranspose,Phi(:,2),'UniformOutput',false);
    end
    PhiX = cell2mat(Phi(:,1).')*kron(speye(size(Phi,1)),X)*cell2mat(Phi(:,2));
else % the superoperator was given as a Choi matrix
    sX = size(X);
    sNX = size(Phi)./sX;
    PhiX = kron(reshape(X,1,prod(sX)),speye(sNX(1)))*reshape(Swap(Phi.',[1,2],[sX(2) sNX(2);sX(1) sNX(1)],1).',sNX(1)*prod(sX),sNX(2));
end
if(~issparse(X))
    PhiX = full(PhiX);
end