%%  CHANNELDISTINGUISHABILITY    Computes the maximum probability of distinguishing two quantum channels
%   This function has two required input arguments:
%     PHI,PSI: linear maps (represented as either Choi matrices or cells of Kraus operators)
%
%   DIST = ChannelDistinguishability(PHI,PSI) is the maximum probability of
%   distinguishing the quantum channels PHI and PSI. Each of PHI and PSI
%   can either be represented by a Choi matrix or by a cell containing
%   Kraus operators of the channel.
%
%   This function has two optional arguments:
%     P (default [1/2, 1/2])
%     DIM (default tries to guess the input and output dimensions of PHI,PSI)
%
%   DIST = ChannelDistinguishability(PHI,PSI,P,DIM) is the maximum
%   probability of distinguishing the quantum channels PHI and PSI, chosen
%   with probability P(1) and P(2), respectively (by default, the channels
%   are chosen uniformly at random). DIM is a 1-by-2 vector containing the
%   input and output dimensions of PHI and PSI (DIM is required if and only
%   if PHI and PSI are provided as Choi matrices and the input and output
%   dimensions are unequal).
%
%   URL: http://www.qetlab.com/ChannelDistinguishability

%   requires: cvx (http://cvxr.com/cvx/), ApplyMap.m, ChoiMatrix.m,
%             ComplementaryMap.m, DiamondNorm.m, DualMap.m, iden.m, IsCP.m,
%             IsHermPreserving.m, IsPSD.m, KrausOperators.m,
%             MaxEntangled.m, opt_args.m, PermuteSystems.m, Swap.m
%
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: November 27, 2014

function dist = ChannelDistinguishability(Phi,Psi,varargin)

% Get the input and output dimensions of PHI and PSI.
[da,db] = superoperator_dims(Psi,0,varargin{2:end});

% set optional argument defaults: p = [1/2,1/2], dim guessed based on input
[p,dim] = opt_args({ [1,1]/2, [da,db] },varargin{:});

% We convert to a Choi matrix to make things easier -- we end up converting
% back to Kraus operators in the DiamondNorm function, but the performance
% impact is negligible, and we have to perform the Kraus decomposition
% again anyway to remove linear dependencies in the Kraus operators, or
% things slow down to a halt in the SDP.
Phi = ChoiMatrix(Phi);
Psi = ChoiMatrix(Psi);

if(length(Phi) ~= length(Psi))
    error('ChannelDistinguishability:DifferentDims','The channels PHI and PSI must have the same dimension input and output spaces as each other.');
elseif(prod(dim) ~= length(Phi))
    error('ChannelDistinguishability:InvalidDim','Could not determine the input and output dimensions of the channels: please provide the DIM argument.');
end

if(abs(sum(p) - 1) > 10*eps || length(p) ~= 2)
    error('ChannelDistinguishability:InvalidP','The vector P must be a probability distribution with 2 entries: its elements must be non-negative and they must sum to 1.');
end

if(max(p) >= 1) % of course we can distinguish 1 object
    dist = 1;
    return
end

% Finally, compute the distinguishability.
dist = DiamondNorm(p(1)*Phi - p(2)*Psi,dim);