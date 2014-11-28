%%  RANDOMSUPEROPERATOR    Generates a random superoperator
%   This function has one required argument:
%     DIM: either a scalar or a 1-by-2 vector specifying the input and
%          output dimensions of the superoperator, in that order
%
%   PHI = RandomSuperoperator(DIM) generates the Choi matrix of a random
%   quantum channel (i.e., a completely positive, trace-preserving linear
%   map) acting on DIM-by-DIM matrices.
%
%   This function has four optional arguments:
%     TP (default 1)
%     UN (default 0)
%     RE (default 0)
%     KR (default prod(DIM))
%
%   PHI = RandomSuperoperator(DIM,TP,UN,RE,KR) generates the Choi matrix of
%   a random superoperator that is trace-preserving if TP=1, is unital if
%   UN=1, has all real entries if RE=1, and has KR or fewer Kraus operators
%   (it will have exactly KR Kraus operators with probability 1).
%
%   URL: http://www.qetlab.com/RandomSuperoperator

%   requires: iden.m, MaxEntangled.m, OperatorSinkhorn.m, opt_args.m,
%             PartialTrace.m, PermuteSystems.m, RandomDensityMatrix.m,
%             RandomStateVector.m, RandomUnitary.m, SchmidtDecomposition.m,
%             Swap.m
%             
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: September 30, 2014

function Phi = RandomSuperoperator(dim,varargin)

% allow the user to enter a single number for dim
if(length(dim) == 1)
    dim = [dim,dim];
end
pd = prod(dim);

% set optional argument defaults: tp=1, un=0, re=0, kr=prod(dim)
[tp,un,re,kr] = opt_args({ 1, 0, 0, pd },varargin{:});

if(tp == 1 && un == 1 && dim(1) ~= dim(2))
    warning('RandomSuperoperator:InvalidDim','There does not exist a unital, trace-preserving map in the case when the input and output dimensions are unequal. The identity matrix will map to a *multiple* of the identity matrix.');
end

% There is a probability 0 chance that the operator Sinkhorn iteration will
% get cranky. Thus we repeatedly try until we *don't* get an error
% (honestly, I'm being slightly overly cautious).
sing_err = 1;
while sing_err
    sing_err = 0;
    try
        % Generate the Choi matrix of a superoperator that is not
        % necessarily trace-preserving or unital. We will enforce those
        % conditions in a moment.
        Phi = RandomDensityMatrix(pd,re,kr,'haar');

        % Now set the appropriate partial traces to the identity.
        if(tp == 1 && un == 0)
            PT = kron(sqrtm(inv(PartialTrace(Phi,2,dim))),eye(dim(2)));
            Phi = PT*Phi*PT;
        elseif(tp == 0 && un == 1)
            PT = kron(eye(dim(1)),sqrtm(inv(PartialTrace(Phi,1,dim))));
            Phi = PT*Phi*PT;
        elseif(tp == 1 && un == 1)
            Phi = OperatorSinkhorn(Phi,dim)*dim(1);
        end
    catch err
        % Operator Sinkhorn error? Try, try again!
        if(strcmpi(err.identifier,'OperatorSinkhorn:LowRank'))
            sing_err = 1;
        else
            rethrow(err);
        end
    end
end