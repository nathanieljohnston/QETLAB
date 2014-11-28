%%  CHOIMAP    Produces the Choi map or one of its generalizations
%   This function has no required arguments.
%
%   C = ChoiMap() is the Choi matrix of the Choi map, which is a positive
%   map on 3-by-3 matrices that is capable of detecting some entanglement
%   that the transpose map is not.
%
%   This function has three optional arguments:
%     A,B,C (default 1,1,0)
%   
%   C = ChoiMap(A,B,C) is the Choi matrix of the positive map defined in
%   [1]. Many of these maps are capable of detecting PPT entanglement.
%
%   URL: http://www.qetlab.com/ChoiMap
%
%   References:
%   [1] S. J. Cho, S.-H. Kye, and S. G. Lee, Linear Algebr. Appl. 171, 213
%   (1992).

%   requires: iden.m, MaxEntangled.m, opt_args.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: August 5, 2013

function C = ChoiMap(varargin)

% set optional argument defaults: a=1, b=1, c=0 (the usual Choi map)
[a,b,c] = opt_args({ 1, 1, 0 },varargin{:});

psi = MaxEntangled(3,0,0);
C = diag([a+1,c,b,b,a+1,c,c,b,a+1]) - psi*psi';