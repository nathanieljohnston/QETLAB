%%  ISKINCOHERENT    Determines whether or not a quantum state is k-incoherent
%   This function has two required arguments:
%     X: a density matrix
%     k: a positive integer, the coherence level to check for
%
%
%   IKI = ISKINCOHERENT(X, K) returns 1 if X is k-incoherent and 0 otherwise.
%
%   This function has no optional arguments.
%
%   URL: http://www.qetlab.com/IskIncoherent
%
%   References:
%   [1] M. Ringbauer, T. R. Bromley, M. Cianciaruso, L. Lami, W. Y. S. Lau, G.
%       Adesso, A. G. White, A. Fedrizzi, and M. Piani. Certification and
%       Quantification of Multilevel Quantum Coherence. Physical Review X 8
%       (2018), no. 4, 041007. https://doi.org/10.1103/PhysRevX.8.041007.
%   [2] N. Johnston, S. Moein, and S. Plosker. The factor width rank of a
%       matrix. Linear Algebra and its Applications 716 (2025), 32–59.
%       https://doi.org/10.1016/j.laa.2025.03.016.

%   requires: has_band_k_ordering.m
%   authors: Benjamin Talbot
%            Luis M. B. Varona (lm.varona@outlook.com)
%   package: QETLAB
%   last updated: July 21, 2025

function iki = IskIncoherent(X, k)
    iki = -1;

    if k <= 0
        error('k must be a positive integer');
    end

    d = size(X,1);

    % Trivial: every quantum state is d-incoherent
    if k == d
        iki = 1;
        return
    end

    if isdiag(X)
        iki = 1;
        return
    end

    % Trivial: only the diagonal quantum states are 1-incoherent
    if k == 1
        iki = 0;
        return
    end

    % [1] Theorem 1
    M = ComparisonMatrix(X);
    if IsPSD(M)
        iki = 1;
        return
    end
    if k == 2
        iki = 0;
        return
    end

    % [1] (8) (Corollary of Theorem 1)
    if trace(X^2) <= 1/(d-1) % and k > 2
        iki = 1;
        return
    end

    % [2] Theorem 1
    if has_band_k_ordering(X, k)
        iki = 1;
        return
    end

    % [1] (7)
    test = (d - k) / (d - 1) * ApplyMap(X, DephasingChannel(d));
    if IsPSD(X - test)
        iki = 1;
        return
    end

    % Hierarchy of the sets of k-incoherence states
    if k >= 2
        iki = IskIncoherent(X, k-1);
    end
    if iki == -1
        iki = KIncohSemidefCheck(X, k);
    end
end

%%  KIncohSemidefCheck    Determines whether or not a quantum state is k-incoherent
%   This function has two required arguments:
%     RHO: a mixed quantum state
%     K: a positive integer
%
%   IKINC = KIncohSemidefCheck(RHO,K) returns 1 if RHO is K-incoherent and 0
%   otherwise. This is checked via semidefinite programming.

%   requires: CVX (http://cvxr.com/cvx/)
%   author: Nathaniel Johnston (nathaniel@njohnston.ca), based on code and
%           conversations with Bartosz Regula
%   last updated: July 21, 2025

function ikinc = KIncohSemidefCheck(rho,k)
    n = size(rho,1);

    Pk = nchoosek(1:n,k);
    s = size(Pk,1);

    cvx_begin sdp quiet
    cvx_precision best;
    variable A(k,k,s) hermitian;
    subject to
        P = zeros(n);

        for j=1:s
            proj = zeros(k,n);

            for i = 1:k
                proj(i,Pk(j,i)) = 1;
            end

            P = P + proj'*A(:,:,j)*proj;

            A(:,:,j) >= 0;
        end

        rho == P;
    cvx_end

    ikinc = 1-min(cvx_optval,1);
end

% Helper function to compute the comparison matrix
% defined as M(i,j) =  |X(i,j)| if i == j, and
%                     -|X(i,j)| if i ~= j
function M = ComparisonMatrix(A)
    dim = size(A);
    M = zeros(dim);

    for i = 1:dim(1)
        for j = 1:dim(2)
            if i == j
                M(i,j) = abs(A(i,j));
            else
                M(i,j) = -abs(A(i,j));
            end
        end
    end
end