%%  HORODECKISTATE    Produces a Horodecki state
%   This function has one required input argument:
%     A: a real parameter in [0,1]
%
%   HORO_STATE = HorodeckiState(A) returns the 3x3 bound entangled
%   Horodecki state described in [1].
%
%   This function has one optional input argument:
%     DIM (default is [3,3], but can be either [3,3] or [2,4])
%
%   HORO_STATE = HorodeckiState(A,DIM) returns the Horodecki state in
%   either (3 \otimes 3)-dimensional space, or (2 \otimes 4)-dimensional
%   space, depending on the dimensions in the 1-by-2 vector DIM.
%
%   The Horodecki state was introduced in [1] which serves as an example in 
%   C^3 \otimes C^3 or C^2 \otimes C^4 of an entangled state that is
%   positive under partial transpose (PPT). The state is PPT for all 
%   a \in [0,1], and separable only for a = 0 or a = 1.
%
%   Note: Refer to [2] (specifically equations (1) and (2)) for more 
%   information on this state and its properties. The 3x3 Horodecki state 
%   is defined explicitly in Section 4.1 of [1] and the 2x4 Horodecki state 
%   is defined explicitly in Section 4.2 of [1].
%
%   References:
%   [1] P. Horodecki. Separability criterion and inseparable mixed states
%       with positive partial transposition. E-print: 
%       arXiv:quant-ph/9703004, 1997.
%       
%   [2] K. Chruscinski. On the symmetry of the seminal Horodecki state.
%       E-print: arXiv:1009.4385 [quant-ph], 2010.
%
%   URL: http://www.qetlab.com/HorodeckiState

%   requires: opt_args.m
%   authors: Vincent Russo (vrusso@uwaterloo.ca)
%            Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB 
%   last updated: December 15, 2014

function horo_state = HorodeckiState( a, varargin )

% set optional argument defaults: dim = [3,3]
[dim] = opt_args({ [3,3] },varargin{:});

if a < 0 || a > 1
    error('HorodeckiState:InvalidA','Argument A must be in the interval [0,1].');
end

if isequal(dim(:),[3;3])
    N_a = 1/(8*a+1);
    b = (1+a)/2;
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
               
elseif isequal(dim(:),[2;4])
    N_a = 1/(7*a+1);
    b = (1+a)/2;
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
    error('HorodeckiState:InvalidDim','DIM must be one of [3,3] or [2,4].');
end