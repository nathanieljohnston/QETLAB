%%  BCSGAMELB    Computes a lower bound on the quantum value of a binary constraint system (BCS) game
%   This function has two required input arguments:
%     D: the local dimension (e.g., D = 2 corresponds to Alice and Bob each
%        having access to a qubit)
%     C: A cell, each of whose elements is a constraint in the BCS. The
%        constraints themselves are specified as 2-by-2-by-2-by-... binary
%        arrays, where the (i,j,k,...)-entry is 1 if and only if setting
%        v1=i, v2=j, v3=k, ... satisfies that constraint.
%
%   BCSLB = BCSGameLB(D,C) is a lower bound on the maximum value that 
%   the specified binary constraint system (BCS) game can take on in
%   quantum mechanical settings where Alice and Bob each have access to
%   D-dimensional quantum systems.
%
%   This function works by starting with a randomly-generated POVM for Bob,
%   and then optimizing Alice's POVM and the shared entangled state. Then
%   Alice's POVM and the entangled state are fixed and Bob's POVM is
%   optimized. And so on, back and forth between Alice and Bob until
%   convergence is reached.
%
%   This function has one optional input argument:
%     VERBOSE (default 1): a flag (either 1 or 0) indicating that the
%     function should or should not print out partial progress as it works. 
%
%   URL: http://www.qetlab.com/BCSGameLB

%   requires: CVX (http://cvxr.com/cvx/), bcs_to_nonlocal.m, opt_args.m,
%             NonlocalGameLB.m
%
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: June 24, 2015

function bcslb = BCSGameLB(d,C,varargin)

    % set optional argument defaults: VERBOSE=1
    [verbose] = opt_args({ 1 },varargin{:});

    % Convert the BCS description of the game to the more general non-local
    % description of the game.
    [p,V] = bcs_to_nonlocal(C);
    
    % Compute a lower bound on the value of this non-local game.
    bcslb = NonlocalGameLB(d,p,V,verbose);
end