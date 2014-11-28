%%  CHESSBOARDSTATE    Produces a chessboard state
%   This function has six required arguments:
%     A,B,C,D,M,N: parameters of the chessboard state, as in [1]
%
%   RHO = ChessboardState(A,B,C,D,M,N) is the chessboard state defined in
%   [1], with S = A*conj(C)/conj(N) and T = A*D/M. If C*M*conj(N) ~=
%   A*B*conj(C) then RHO is entangled. If each of A,B,C,D,M,N are real then
%   RHO has positive partial transpose, and is hence bound entangled.
%
%   This function has two optional arguments:
%     S (default A*conj(C)/conj(N))
%     T (default A*D/M)
%   
%   RHO = ChessboardState(A,B,C,D,M,N,S,T) is the chessboard state defined
%   in [1]. Note that, for certain choices of S and T, this state will not
%   have positive partial transpose, and thus may not be bound entangled --
%   a warning will be produced in these cases.
%
%   URL: http://www.qetlab.com/ChessboardState
%
%   References:
%   [1] D. Bruss and A. Peres. Construction of quantum states with bound
%       entanglement. Phys. Rev. A, 61:30301(R), 2000.

%   requires: IsPPT.m, IsPSD.m, opt_args.m, PartialTranspose.m,
%             PermuteSystems.m
%
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: March 13, 2013

function rho = ChessboardState(a,b,c,d,m,n,varargin)

% set optional argument defaults: s = ac*/n*, t = ad/m
[s,t] = opt_args({ a*conj(c)/conj(n), a*d/m },varargin{:});

v1 = [m,0,s,0,n,0,0,0,0];
v2 = [0,a,0,b,0,c,0,0,0];
v3 = [conj(n),0,0,0,-conj(m),0,t,0,0];
v4 = [0,conj(b),0,-conj(a),0,0,0,d,0];

rho = v1'*v1 + v2'*v2 + v3'*v3 + v4'*v4;
rho = rho/trace(rho);

if(~IsPPT(rho))
    warning('ChessboardState:NotPPT','The specified chessboard state does not have positive partial transpose.');
end