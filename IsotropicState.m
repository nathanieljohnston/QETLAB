%%  ISOTROPICSTATE    Produces an isotropic state
%   This function has two required arguments:
%     DIM: the local dimension
%     ALPHA: the parameter of the isotropic state
%
%   RHO = IsotropicState(DIM,ALPHA) is the isotropic state with parameter
%   ALPHA acting on (DIM*DIM)-dimensional space. More specifically, RHO is
%   the density operator defined by (1-ALPHA)*I/DIM^2 + ALPHA*E, where I is
%   the identity operator and E is the projection onto the standard
%   maximally-entangled pure state on two copies of DIM-dimensional space.
%
%   URL: http://www.qetlab.com/IsotropicState

%   requires: iden.m, MaxEntangled.m, opt_args.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: September 22, 2014

function rho = IsotropicState(dim,alpha)

% compute the isotropic state
psi = MaxEntangled(dim,1,0);
rho = (1-alpha)*speye(dim^2)/dim^2 + alpha*psi*psi'/dim;