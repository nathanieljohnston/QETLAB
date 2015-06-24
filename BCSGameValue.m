%%  BCSGAMEVALUE    Computes the maximum value of a binary constraint system (BCS) game
%   This function has one required input argument:
%     C: A cell, each of whose elements is a constraint in the BCS. The
%        constraints themselves are specified as 2-by-2-by-2-by-... binary
%        arrays, where the (i,j,k,...)-entry is 1 if and only if setting
%        v1=i, v2=j, v3=k, ... satisfies that constraint.
%
%   BCSVAL = BCSGameValue(C) is the maximum value that the specified 
%   BCS game can take on in classical mechanics. For the maximum quantum
%   or no-signalling value, see the optional arguments described below.
%
%   This function has two optional input arguments:
%     MTYPE (default 'classical'): one of 'classical', 'quantum', or
%       'nosignal', indicating what type of BCS game value should be
%       computed. IMPORTANT NOTE: if MTYPE='quantum' then only an upper
%       bound on the nonlocal game value is computed, not its exact value
%       (see the argument K below).
%     K (default 1): if MYTPE='quantum', then K is a non-negative integer
%       or string indicating what level of the NPA hierarchy to use to
%       bound the nonlocal game (higher values give better bounds, but
%       require more computation time). See the NPAHierarchy function for
%       details.
%
%   BCSVAL = BCSGameValue(C,MTYPE,K) is the maximum value that the
%   specified nonlocal game can take on in the setting (classical, quantum,
%   or no-signalling) specified by MTYPE.
%
%   URL: http://www.qetlab.com/BCSGameValue

%   requires: CVX (http://cvxr.com/cvx/), bcs_to_nonlocal.m,
%             NonlocalGameValue.m, opt_args.m
%
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: June 24, 2015

function bcsval = BCSGameValue(C,varargin)

    % set optional argument defaults: MTYPE='classical', K=1
    [mtype,k] = opt_args({ 'classical', 1 },varargin{:});

    % Convert the BCS description of the game to the more general non-local
    % description of the game.
    [p,V] = bcs_to_nonlocal(C);
    
    % Compute the value of this non-local game.
    bcsval = NonlocalGameValue(p,V,mtype,k);
end