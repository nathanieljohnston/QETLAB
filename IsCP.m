%%  ISCP    Determines whether or not a superoperator is completely positive
%   This function has one required argument:
%     PHI: a superoperator
%
%   CP = IsCP(PHI) is either 1 or 0, indicating that PHI is or is not
%   completely positive (within reasonable numerical error).
%
%   This function has one optional input argument:
%     TOL (default eps^(3/4))
%
%   CP = IsCP(PHI,TOL) determines whether or not PHI is completely positive
%   within the numerical tolerance specified by TOL.
%
%   URL: http://www.qetlab.com/IsCP

%   requires: ApplyMap.m, ChoiMatrix.m, iden.m, IsHermPreserving.m,
%             IsPSD.m, MaxEntangled.m, opt_args.m, PermuteSystems.m,
%             sporth.m, superoperator_dims.m
%
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: January 4, 2013

function cp = IsCP(Phi,varargin)

if(iscell(Phi) && size(Phi,2) == 1)
    cp = 1;
    return
end

% Use Choi's theorem to determine whether or not Phi is CP.
C = ChoiMatrix(Phi);
cp = (IsHermPreserving(C,varargin{:}) && IsPSD(C,varargin{:}));