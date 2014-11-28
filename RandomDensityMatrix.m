%%  RANDOMDENSITYMATRIX    Generates a random density matrix
%   This function has one required argument:
%     DIM: the number of rows (and columns) of the density matrix
%
%   RHO = RandomDensityMatrix(DIM) generates a random DIM-by-DIM density
%   matrix, distributed according to the Hilbert-Schmidt measure.
%
%   This function has three optional arguments:
%     RE (default 0)
%     K (default DIM)
%     DIST (default 'haar')
%
%   RHO = RandomDensityMatrix(DIM,RE,K,DIST) generates a random density
%   matrix of rank <= K, distributed according to the distribution DIST. If
%   RE = 1 then all of its entries will be real. DIST must be one of:
%     'haar' or 'hs' (default) - Generate a larger  pure state according to
%                      Haar measure and trace out the extra dimensions.
%                      Sometimes called the Hilbert-Schmidt measure when
%                      K = DIM.
%     'bures'        - the Bures measure
%
%   URL: http://www.qetlab.com/RandomDensityMatrix

%   requires: opt_args.m, RandomUnitary.m
%
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: October 10, 2014

function rho = RandomDensityMatrix(dim,varargin)

% set optional argument defaults: re=0, k=dim, dist='hs'
[re,k,dist] = opt_args({ 0, dim, 'haar' },varargin{:});

% Haar/Hilbert-Schmidt measure
gin = randn(dim,k);
if(~re)
    gin = gin + 1i*randn(dim,k);
end
if(strcmpi(dist,'bures')) % Bures measure
    gin = (RandomUnitary(dim,re) + eye(dim))*gin;
end

rho = gin*gin';
rho = rho/trace(rho);