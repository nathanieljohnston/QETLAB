%%	GISINSTATE		Produces a Gisin state
%	This function has two required input argument:
%       LAMBDA: Parameter a \in [0,1] is responsible for calculating the 
%          resultant matrix.
%       THETA: 
%
%   GISIN_STATE = GisinState(LAMBDA, THETA) returns the Gisin state
%   described in [1].
%
%   The Gisin states are a mixture of the entangled state rho_theta and the 
%   separable states rho_uu and rho_dd. 
%   
%	References:
%   [1] Hidden quantum nonlocality revealed by local filters. N. Gisin.
%       (http://dx.doi.org/10.1016/S0375-9601(96)80001-6). 1996.
%
%	URL: http://www.qetlab.com/GisinState

%	requires: Nothing
%
% 	author: Vincent Russo (vrusso@uwaterloo.ca)
%           Nathaniel Johnston (nathaniel@njohnston.ca)
%	package: QETLAB 
%	last updated: November 27, 2014

function gisin_state = GisinState( lambda, theta )

if lambda < 0 || lambda > 1
    error('GisinState:ImproperVal','Lambda must be 0 <= lambda <= 1.');
end

rho_theta = [ 0     0                0                0;
              0   sin(theta)^2      -1/2*sin(2*theta) 0;
              0   -1/2*sin(2*theta) cos(theta)^2      0;
              0     0                0                0 ];
          
rho_uu = [ 1 0 0 0;
           0 0 0 0;
           0 0 0 0;
           0 0 0 0 ];
       
rho_dd = [ 0 0 0 0;
           0 0 0 0;
           0 0 0 0;
           0 0 0 1 ];

gisin_state = lambda * rho_theta + 1/2*(1-lambda)*(rho_uu + rho_dd);

end