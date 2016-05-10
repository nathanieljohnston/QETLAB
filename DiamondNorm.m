%%  DIAMONDNORM    Computes the diamond norm of a superoperator
%   This function has one required input argument:
%     PHI: a superoperator
%
%   DN = DiamondNorm(PHI) is the diamond norm of the superoperator PHI. PHI
%   should be specified as either a cell with one or two columns of Kraus
%   operators, or as a Choi matrix (see online QETLAB tutorial for details
%   about specifying superoperators). If PHI is provided as a Choi matrix
%   with unequal input and output dimensions, a second argument specifying
%   the dimensions should also be provided (see below).
%
%   This function has one optional input argument:
%     DIM (default has both subsystems of equal dimension)
%
%   DN = DiamondNorm(PHI,DIM) is the diamond norm of PHI, as above, where
%   DIM is a 1-by-2 vector containing the input and output dimensions of
%   PHI, in that order (equivalently, these are the dimensions of the first
%   and second subsystems of the Choi matrix PHI, in that order).
%
%   URL: http://www.qetlab.com/DiamondNorm

%   requires: ApplyMap.m, ComplementaryMap.m, ChoiMatrix.m, CVX
%             (http://cvxr.com/cvx/), DualMap.m, iden.m, IsCP.m,
%             IsHermPreserving.m, IsPSD.m, KrausOperators.m,
%             MaxEntangled.m, opt_args.m, PermuteSystems.m, Swap.m
%
%   author: Nathaniel Johnston (nathaniel@njohnston.ca), based on an
%           algorithm by John Watrous
%   package: QETLAB
%   last updated: April 26, 2016

function dn = DiamondNorm(Phi,varargin)

% If PHI is a CVX variable, assume it is already a Choi matrix and get its
% dimensions in a CVX-safe way.
if(isa(Phi,'cvx'))
    len_phi = length(Phi);
    sqrt_len = round(sqrt(len_phi));

    % Set optional argument defaults: dim=sqrt_len
    [dim] = opt_args({ [sqrt_len,sqrt_len] },varargin{:});
    
% If PHI is not a CVX variable, we can get its dimensions a bit more
% robustly and allow Kraus operator input. We can also implement some
% speed-ups like computing the diamond norm of a completely positive map
% exactly.
else
    Phi = KrausOperators(Phi,varargin{:});
    [dim] = opt_args({ [size(Phi{1,1},2),size(Phi{1,1},1)] },varargin{:});
    len_phi = dim(1)*dim(2);

    two_cols = (size(Phi,2) > 1);
    if(two_cols && ((size(Phi{1,2},2) ~= dim(1) || size(Phi{1,2},1) ~= dim(2))))
        error('DiamondNorm:InvalidDim','The input and output spaces of PHI must both be square.');
    end

    % The diamond norm of a CP map is trivial to compute.
    if(~two_cols || IsCP(Phi))
        dn = norm(ApplyMap(speye(dim(2)),DualMap(Phi)));
        return
    end
    
    Phi = ChoiMatrix(Phi);
end

% This SDP has two advantages over the other SDP for the diamond norm:
% it works in a CVX-safe way that allows this function to be used within
% other CVX optimization problems, and it has better numerical accuracy.
% The downside is that it is slightly slower for channels with few Kraus
% operators.
cvx_begin sdp quiet
    cvx_precision best;
    variable Y0(len_phi,len_phi) hermitian
    variable Y1(len_phi,len_phi) hermitian
    minimize norm(PartialTrace(Y0,2,dim)) + norm(PartialTrace(Y1,2,dim))
    subject to
        cons = [Y0,-Phi;-Phi',Y1];
        cons + cons' >= 0; % avoid some numerical problems: CVX often thinks things aren't symmetric without this
        Y0 >= 0;
        Y1 >= 0;
cvx_end
    
dn = cvx_optval/2;