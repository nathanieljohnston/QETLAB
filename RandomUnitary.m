%%  RANDOMUNITARY    Generates a random unitary or orthogonal matrix
%   This function has one required argument:
%     DIM: the number of rows (and columns) of the unitary matrix
%
%   U = RandomUnitary(DIM) generates a random DIM-by-DIM unitary matrix,
%   uniformly distributed according to Haar measure.
%
%   This function has one optional argument:
%     RE (default 0)
%
%   U = RandomUnitary(DIM,RE) generates a random unitary matrix (if
%   RE=0) or a random real orthogonal matrix (if RE=1), uniformly
%   distributed according to the Haar measure.
%
%   URL: http://www.qetlab.com/RandomUnitary

%   requires: opt_args.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: September 30, 2014

function U = RandomUnitary(dim,varargin)

% set optional argument defaults: re=0
[re] = opt_args({ 0 },varargin{:});

% construct the Ginibre ensemble
gin = randn(dim);
if(~re)
    gin = gin + 1i*randn(dim);
end

% QR decomposition of the Ginibre ensemble
[Q,R] = qr(gin);

% compute U from the QR decomposition
R = sign(diag(R));
R(R==0) = 1; % protect against potentially zero diagonal entries
U = bsxfun(@times,Q,R.'); % much faster than the naive U = Q*diag(R)