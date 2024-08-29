%%  COHERENCERANK    Calculates the coherence rank of a pure state
%   This function has one required argument:
%     v: the vector for which to calculate the coherence rank, assumed to be
%        represented in the standard basis
%
%   cohRank = CoherenceRank(v) calculates the coherence rank of the pure
%   state v with respect to the standard basis. The coherence rank is defined
%   as the number of non-zero entries in v when v is expressed in a given basis.
%
%   This function has two optional arguments:
%     tol (default 1e-10)
%     basis (default standard basis)
%   
%   cohRank = CoherenceRank(X,TOL,BASIS,@INNERPRODUCT) calculates the coherence 
%   rank of the pure state v with respect to the basis BASIS, using a tolerance
%   of TOL. BASIS should be a unitary matrix whose columns are the basis vectors.
%
%   URL: http://www.qetlab.com/CoherenceRank
%   
%   References:
%   [1] Ringbauer, Martin and Bromley, Thomas R. and Cianciaruso, Marco and
%       Lami, Ludovico and Lau, W. Y. Sarah and Adesso, Gerardo and White,
%       Andrew G. and Fedrizzi, Alessandro and Piani, Marco. Certification
%       and Quantification of Multilevel Quantum Coherence. American Physical 
%       Society, 10.1103/PhysRevX.8.041007, 2018.

%   requires: opt_args.m
%   author: Benjamin Talbot
%   package: QETLAB
%   last updated: August 25, 2024

function cohRank = CoherenceRank(v, varargin)
    [basis, tol] = opt_args({speye(size(v)), 1e-10}, varargin{:});
    if (isrow(v))
        v = v';
    end
    if length(varargin) > 1
        v = basisCoordinates();
    end

    cohRank = 0;
    for i = 1:size(v)
        if abs(v(i)) <= tol
            cohRank = cohRank + 1;
        end
    end

    % Compute the basis coordinates of the input
    % vector with respect to a given basis
    function basisCoords = basisCoordinates()
        basisCoords = basis \ v;
    end
end