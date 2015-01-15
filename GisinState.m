%%  GISINSTATE    Produces a Gisin state
%   This function has two required input argument:
%     LAMBDA: A real parameter in [0,1].
%     THETA: A real parameter.
%
%   GISIN_STATE = GisinState(LAMBDA,THETA) returns the Gisin state
%   described in [1].
%
%   The Gisin states are a mixture of the entangled state rho_theta and the 
%   separable states rho_uu and rho_dd. 
%   
%   References:
%   [1] N. Gisin. Hidden quantum nonlocality revealed by local filters.
%       (http://dx.doi.org/10.1016/S0375-9601(96)80001-6). 1996.
%
%   URL: http://www.qetlab.com/GisinState

%   requires: nothing
%   author: Vincent Russo (vrusso@uwaterloo.ca)
%   package: QETLAB 
%   last updated: January 14, 2015

function gisin_state = GisinState( lambda, theta )

if lambda < 0 || lambda > 1
    error('GisinState:ImproperVal','LAMBDA must satisfy 0 <= LAMBDA <= 1.');
end

rho_theta = [ 0     0                0                0;
              0   sin(theta)^2      -sin(2*theta)/2 0;
              0   -sin(2*theta)/2   cos(theta)^2      0;
              0     0                0                0 ];
          
rho_uu_dd = [ 1 0 0 0;
              0 0 0 0;
              0 0 0 0;
              0 0 0 1 ];

gisin_state = lambda*rho_theta + (1-lambda)*rho_uu_dd/2;

end