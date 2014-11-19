%%	HORODECKISTATE		Produces a Horodecki state
%	This function has one required input argument:
%       A: Parameter a \in [0,1] is responsible for calculating the 
%          resultant matrix.
%
%   HORO_STATE = HorodeckiState(A) returns the 3x3 Horodecki state is no
%   dimension is specified.
%
%   This function has one optional input argument:
%      DIM: (default is 3x3, but can be either 2x3 or 3x3)
%
%   HORO_STATE = HorodeckiState(A,DIM) return
%
%   The Horodecki state was introduced in [1] which serves as an example in 
%   C^3 \otimes C^3 or C^2 \otimes C^4 which represents an entangled state 
%   positive under partial transpose (PPT). The state is PPT for all 
%   a \in [0,1], and separable for a = 0 or a = 1.
%
%   Note: Refer to [2] (specifically equations (1) and (2)) for more 
%   information on this state and its properties. The 3x3 Horodecki state 
%   is defined explicitly in Section 4.1 in [1] and the 2x4 Horodecki state 
%   is definite explicitly in Section 4.2 in [1].
%
%	References:
%   [1] P. Horodecki. Separability criterion and inseparable mixed states
%       with positive partial transposition. E-print: 
%       arXiv:quant-ph/9703004, 1997.
%       
%   [2] K. Chruscinski. On the symmetry of the seminal Horodecki state.
%       E-print: arXiv:1009.4385 [quant-ph], 2010.
%
%	URL: http://www.qetlab.com/HorodeckiState

%	requires: Nothing
%
% 	author:
%	package: 
%	version: 
%	last updated:

function horo_state = HorodeckiState( a, varargin )

if a < 0 || a > 1
    error('HorodeckiState:ImproperDim','Argument a must be  in[0,1].');
end

if nargin >= 2
    dim = varargin{1};
else
    dim = '3x3';
end

if strcmp(dim,'3x3')
    N_a = 1/(8*a+1);
    b = 1+a/2;
    c = sqrt(1-a^2)/2;

    horo_state = N_a * [ a 0 0 0 a 0 0 0 a;
                         0 a 0 0 0 0 0 0 0;
                         0 0 a 0 0 0 0 0 0;
                         0 0 0 a 0 0 0 0 0;
                         a 0 0 0 a 0 0 0 a;
                         0 0 0 0 0 a 0 0 0;
                         0 0 0 0 0 0 b 0 c;
                         0 0 0 0 0 0 0 a 0;
                         a 0 0 0 a 0 c 0 b ];
               
elseif strcmp(dim, '2x4')
    N_a = 1/(7*a+1);
    b = 1+a/2;
    c = sqrt(1-a^2)/2;
    
    horo_state = N_a * [ a 0 0 0 0 a 0 0;
                         0 a 0 0 0 0 a 0;
                         0 0 a 0 0 0 0 a;
                         0 0 0 a 0 0 0 0;
                         0 0 0 0 b 0 0 c;
                         a 0 0 0 0 a 0 0;
                         0 a 0 0 0 0 a 0;
                         0 0 a 0 c 0 0 b ];
else
    error('HorodeckiState:ImproperState','Dim must be either 2x4 or 3x3.');
end
        
end