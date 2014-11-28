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
%   last updated: November 21, 2014

function dn = DiamondNorm(Phi,varargin)

% If PHI is a CVX variable, compute the diamond norm in a CVX-safe way that
% allows this function to be used within other CVX optimization problems.
if(isa(Phi,'cvx'))
    len_phi = length(Phi);
    sqrt_len = round(sqrt(len_phi));

    % Set optional argument defaults: dim=sqrt_len
    [dim] = opt_args({ [sqrt_len,sqrt_len] },varargin{:});

    cvx_begin sdp quiet
        cvx_precision default;
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

% If PHI is not a CVX variable, compute the diamond norm via a slightly
% different SDP that will be faster if PHI has few Kraus operators (but
% will require about the same amount of computation time if PHI has a full
% set of Kraus operators). Also implement other non-CVX-save speedups, like
% quickly computing the diamond norm of a CP map.
else
    % Always call the KrausOperators function, even if Phi is already a set of
    % Kraus operators, to make sure we get a minimal set -- this is important
    % for speed reasons.
    Phi = KrausOperators(Phi,varargin{:});

    da = size(Phi{1,1},2);
    db = size(Phi{1,1},1);
    nk = size(Phi,1);
    two_cols = (size(Phi,2) > 1);
    if(two_cols && ((size(Phi{1,2},2) ~= da || size(Phi{1,2},1) ~= db)))
        error('DiamondNorm:InvalidDim','The input and output spaces of PHI must both be square.');
    end

    % The diamond norm of a CP map is trivial to compute.
    if(~two_cols || IsCP(Phi))
        dn = norm(ApplyMap(speye(db),DualMap(Phi)));
        return
    end

    % Otherwise, compute the diamond norm via semidefinite programming.
    cvx_begin sdp quiet
        cvx_precision default;
        variable rho0(da,da) hermitian
        variable rho1(da,da) hermitian
        variable X(nk,nk) complex
        maximize trace(X) + trace(X')
        subject to
            cons = [ApplyMap(rho0,ComplementaryMap(Phi(:,1))),X;X',ApplyMap(rho1,ComplementaryMap(Phi(:,2)))];
            cons + cons' >= 0; % avoid some numerical problems: CVX often thinks things aren't symmetric without this
            trace(rho0) <= 1/2;
            trace(rho1) <= 1/2;
            rho0 >= 0;
            rho1 >= 0;
    cvx_end

    dn = real(cvx_optval);
end