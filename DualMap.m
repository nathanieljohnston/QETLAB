%%  DUALMAP    Computes the dual of a superoperator in the Hilbert-Schmidt inner product
%   This function has one required argument:
%     PHI: a superoperator
%
%   PHID = DualMap(PHI) is the dual map (in the Hilbert-Schmidt inner
%   product) of the superoperator PHI. If PHI is provided as a cell of
%   Kraus operators, then so is PHID. If PHI is provided as a Choi matrix,
%   then so is PHID. If PHI is provided as a Choi matrix the input space is
%   not of dimension equal to that of the output space, a second argument
%   DIM must be provided (see below).
%
%   This function has one optional input argument:
%     DIM (default has input and output of equal dimension)
%
%   PHID = DualMap(PHI,DIM) is the same as above, where DIM is a
%   1-by-2 vector containing the input and output dimensions of PHI, in
%   that order (equivalently, these are the dimensions of the first and
%   second subsystems of the Choi matrix PHI, in that order). If the input
%   or output space is not square, then DIM's first row should contain the
%   input and output row dimensions, and its second row should contain its
%   input and output column dimensions. DIM is required if and only if PHI
%   has unequal input and output dimensions and is provided as a Choi
%   matrix.
%
%   URL: http://www.qetlab.com/DualMap

%   requires: opt_args.m, PermuteSystems.m, sporth.m, superoperator_dims.m,
%             Swap.m
%
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: November 24, 2014

function PhiD = DualMap(Phi,varargin)

if(iscell(Phi)) % Phi is provided as a set of Kraus operators
    PhiD = cellfun(@ctranspose,Phi,'UniformOutput',false);
else % Phi is provided as a Choi matrix
    % Get the dimensions of PHI.
    [da,db] = superoperator_dims(Phi,1,varargin{:});
    
    % Compute the dual map.
    PhiD = Swap(conj(Phi),[1,2],[da.',db.']);
end