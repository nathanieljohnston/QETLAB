%%	BREUERSTATE		Produces a Breuer State of dimension >= 4
%	This function has two required arguments:
%		DIM: the local dimension (for dim >= 4 and even) 
%		LAMBDA: describes the weight of the singlet component
%
%	BREUER_STATE = BreuerState(DIM, LAMBDA) Gives a Breuer bound entangled 
%                  state for two qubits of dimension DIM, with the LAMBDA 
%                  parameter describing the weight of the singlet 
%                  component. 
%	References:
%   [1] H-P. Breuer. Optimal entanglement criterion for mixed quantum 
%       states. E-print: arXiv:quant-ph/0605036, 2006.
%
%	URL: http://www.qetlab.com/BreuerState

%	requires: SymmetricProjection.m
%
% 	author:
%	package: 
%	version: 
%	last updated: 

function breuer_state = BreuerState(dim, lambda)

% SU(2) generators for dxd matrices
N = dim-1;
a1da1 = diag(0:N);
a2da2 = diag(N-(0:N));
a1da2 = zeros(N+1,N+1);
for k = 0:N-1
    a1da2(k+1+1,k+1) = sqrt((k+1)*(N-k));   
end
a2da1 = a1da2';

% Schwinger's construction
jz = (a1da1-a2da2)/2;
jx = (a1da2+a2da1)/2;
jy = -1i*(a1da2-a2da1)/2;

% Collective operators
ee = eye(dim);
Jx = kron(jx,ee) + kron(ee,jx);
Jy = kron(jy,ee) + kron(ee,jy);
Jz = kron(jz,ee) + kron(ee,jz);

% A singlet is the ground state of the following Hamiltonian
H = Jx^2 + Jy^2 + Jz^2;
[v,D] = eig(H);
[~,index] = min(diag(real(D)));
phis = v(:,index);

% Mix it with the normalized projector to the symmetric subspace
breuer_state = lambda * (phis*phis') + (1-lambda) * ... 
               SymmetricProjection(dim)/trace(SymmetricProjection(dim));

end