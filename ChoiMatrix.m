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

%   requires: ApplyMap.m, iden.m, MaxEntangled.m, opt_args.m, PartialMap.m,
%             PermuteSystems.m, sporth.m, superoperator_dims.m
%
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: November 24, 2014

function C = ChoiMatrix(Phi,varargin)

if(~iscell(Phi)) % PHI is already a Choi matrix
    C = Phi;
    return
end

% Get the dimensions of PHI.
da = superoperator_dims(Phi); % we only need the input space dimensions

% Set optional argument defaults: sys=2
[sys] = opt_args({ 2 },varargin{:});

% Now create the Choi matrix: apply the map to half of the (unnormalized)
% maximally-entangled states.
C = PartialMap(MaxEntangled(da(1),1,0)*MaxEntangled(da(2),1,0).',Phi,sys,[da(1),da(1);da(2),da(2)]);
if(~issparse(Phi{1,1}))
    C = full(C);
end