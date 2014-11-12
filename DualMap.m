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

%   requires: opt_args.m, PermuteSystems.m, Swap.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   version: 0.50
%   last updated: November 12, 2014

function PhiD = DualMap(Phi,varargin)

if(iscell(Phi)) % Phi is provided as a set of Kraus operators
    PhiD = cellfun(@ctranspose,Phi,'UniformOutput',false);
else % Phi is provided as a Choi matrix
    dX = size(Phi);
    sdX = round(sqrt(dX));

    % Set optional argument defaults: dim=sqrt(length(Phi))
    [dim] = opt_args({ [sdX(1) sdX(1);sdX(2) sdX(2)] },varargin{:});

    % allow the user to enter a single number for dim
    if(length(dim) == 1)
        dim = [dim,dX(1)/dim];
        if abs(dim(2) - round(dim(2))) >= 2*dX(1)*eps
            error('DualMap:InvalidDim','If DIM is a scalar, PHI must be square and DIM must evenly divide length(PHI); please provide the DIM array containing the dimensions of the subsystems.');
        end
        dim(2) = round(dim(2));
    end
    
    % allow the user to enter a vector for dim if X is square
    if(min(size(dim)) == 1)
        dim = dim(:).'; % force dim to be a row vector
        dim = [dim;dim];
    end
    
    PhiD = Swap(conj(Phi),[1,2],dim);
end