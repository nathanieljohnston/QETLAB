%% CHOIMATRIX     Computes the Choi matrix of a superoperator
%   This function has one required argument:
%     PHI: a superoperator
%
%   C = ChoiMatrix(PHI) is the Choi matrix of the superoperator PHI. The
%   default convention used within this function (and throughout QETLAB) is
%   that the Choi matrix is the result of applying the map PHI to the
%   *second* subsystem of the standard maximally entangled (unnormalized)
%   state. PHI should be provided either as a Choi matrix, or as a cell
%   with either 1 or 2 columns whose entries are its Kraus operators (see
%   full QETLAB documentation for details).
%
%   This function has one optional argument:
%     SYS (default 2)
%
%   C = ChoiMatrix(PHI,SYS) is the Choi matrix of PHI, with the convention
%   that the map PHI is applied to the SYS subsystem of the standard
%   maximally entangled state. SYS must be either 1 or 2.
%
%   URL: http://www.qetlab.com/ChoiMatrix

%   requires: ApplyMap.m, iden.m, MaxEntangled.m, opt_args.m,
%             PermuteSystems.m
%
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   version: 0.50
%   last updated: January 21, 2013

function C = ChoiMatrix(Phi,varargin)

if(~iscell(Phi)) % is already a Choi matrix
    C = Phi;
    return
end

% set optional argument defaults: sys=2
[sys] = opt_args({ 2 },varargin{:});

sPhi = size(Phi);
n1 = size(Phi{1,1},2);
psi1 = MaxEntangled(n1,1,0);
if(sPhi(2) == 1 || (sPhi(1) == 1 && sPhi(2) > 2)) % map is CP
    n2 = n1;
    psi2 = psi1;
else
    n2 = size(Phi{1,2},2);
    psi2 = MaxEntangled(n2,1,0);
end

C = PartialMap(psi1*psi2.',Phi,sys,[n1,n1;n2,n2]);
if(~issparse(Phi{1,1}))
    C = full(C);
end