%%  FIDELITY    Computes the fidelity of two density matrices
%   This function has two required input arguments:
%     RHO,SIGMA: density matrices
%
%   FID = Fidelity(RHO,SIGMA) is the fidelity between the two density
%   matrices RHO and SIGMA, defined by ||sqrt(RHO)*sqrt(SIGMA)||_1, where
%   ||.||_1 denotes the trace norm. FID is a value between 0 and 1, with 0
%   corresponding to matrices RHO and SIGMA with orthogonal support, and 1
%   corresponding to the case RHO = SIGMA.
%
%   URL: http://www.qetlab.com/Fidelity

%   requires: nothing (only uses CVX if called as an input to CVX)
%   authors: Vincent Russo (vrusso@uwaterloo.ca)
%            Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: November 21, 2014

function fid = Fidelity(rho,sigma)

sz_rho = size(rho);

% Do some error checking.
if(~all(sz_rho == size(sigma)))
    error('Fidelity:InvalidDims','RHO and SIGMA must be matrices of the same size.');
elseif(sz_rho(1) ~= sz_rho(2))
    error('Fidelity:InvalidDims','RHO and SIGMA must be square.');
end

% If rho or sigma is a CVX variable then compute fidelity via semidefinite
% programming, so that this function can be used in the objective function
% or constraints of other CVX optimization problems.
if(isa(rho,'cvx') || isa(sigma,'cvx'))
    cvx_begin sdp quiet
        variable X(sz_rho(1),sz_rho(1)) complex;
        maximize trace(X) + trace(X');
        subject to
            cons = [rho,X;X',sigma];
            cons + cons' >= 0; % avoid some numerical problems: CVX often thinks things aren't symmetric without this
    cvx_end
    
    fid = cvx_optval/2;

% If rho and sigma are *not* CVX variables, compute fidelity normally,
% since this is much faster.
else
    [sq_rho,res] = sqrtm(rho); % need "res" parameter to suppress MATLAB singularity warning
    [sq_fid,res] = sqrtm(sq_rho*sigma*sq_rho);
    fid = real(trace(sq_fid)); % finally, compute the fidelity
end