%%  MATSUMOTOFIDELITY    Computes the Matsumoto fidelity of two density matrices
%   This function has two required input arguments:
%     RHO,SIGMA: density matrices
%
%   FID = MatsumotoFidelity(RHO,SIGMA) is the Matsumoto fidelity between
%   the two density matrices RHO and SIGMA, defined by tr(RHO#SIGMA), where
%   RHO#SIGMA is the matrix geometric mean. FID is a value between 0 and 1,
%   with 0 corresponding to matrices RHO and SIGMA with 0-intersection
%   ranges, and 1 corresponding to the case RHO = SIGMA.
%
%   URL: http://www.qetlab.com/MatsumotoFidelity

%   requires: nothing (only uses CVX if called as an input to CVX)
%   authors: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: June 23, 2020

function fid = MatsumotoFidelity(rho,sigma)

sz_rho = size(rho);

% Do some error checking.
if(~all(sz_rho == size(sigma)))
    error('MatsumotoFidelity:InvalidDims','RHO and SIGMA must be matrices of the same size.');
elseif(sz_rho(1) ~= sz_rho(2))
    error('MatsumotoFidelity:InvalidDims','RHO and SIGMA must be square.');
end

% If rho or sigma is a CVX variable then compute fidelity via semidefinite
% programming, so that this function can be used in the objective function
% or constraints of other CVX optimization problems.
if(isa(rho,'cvx') || isa(sigma,'cvx'))
    cvx_begin sdp quiet
        cvx_precision best
        variable X(sz_rho(1),sz_rho(1)) hermitian;
        maximize trace(X);
        subject to
            cons = [rho,X;X,sigma];
            cons + cons' >= 0; % avoid some numerical problems: CVX often thinks things aren't hermitian without this
    cvx_end
    
    fid = cvx_optval;

% If rho and sigma are *not* CVX variables, compute fidelity normally,
% since this is much faster.
else
    if(abs(det(sigma)) > abs(det(rho))) % for numerical stability, invert the density matrix with larger determinant
        temp_rho = rho;
        rho = sigma;
        sigma = temp_rho;
    end
    rho = rho + 10^(-8)*eye(sz_rho(1)); % need rho to be invertilble -- this is a ham-fisted way of avoiding problems
    [sq_rho,res] = sqrtm(rho); % need "res" parameter to suppress MATLAB singularity warning
    isq_rho = inv(sq_rho);
    isq_rho = (isq_rho+isq_rho')/2;
    [sq_fid,res] = sqrtm(isq_rho*sigma*isq_rho);
    fid = real(trace(sq_rho*sq_fid*sq_rho)); % finally, compute the fidelity
end