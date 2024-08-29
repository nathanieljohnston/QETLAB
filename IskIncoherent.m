%%  ISKINCOHERENT    Determines whether or not a quantum state is k-incoherent
%   This function has two required arguments:
%     X: a density matrix
%     k: a positive integer, the coherence level to check for
%
%   IKI = ISKINCOHERENT(X, K) returns 1 if X is k-incoherent, 0 if X is not
%   k-incoherent, and -1 if the function cannot determine if X is k-incoherent
%   or not.
%   
%   URL: http://www.qetlab.com/IskIncoherent
%   
%   References:
%   [1] Ringbauer, Martin and Bromley, Thomas R. and Cianciaruso, Marco and Lami,
%       Ludovico and Lau, W. Y. Sarah and Adesso, Gerardo and White, Andrew G. 
%       and Fedrizzi, Alessandro and Piani, Marco. Certification and Quantification 
%       of Multilevel Quantum Coherence. Phys. Rev. X, 8.041007, 2018.

%   author: Benjamin Talbot
%   package: QETLAB
%   last updated: August 26, 2024


% DOUBLE CHECK FOR THE COHERENCE NUMBER FROM THIS PAPER, MIGHT THROW OFF SOME OF MY ASSUMPTIONS


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
    % Trivial: only the diagonal quantum states are 1-incoherent
    elseif k == 1
        iki = 0;
        return
    end

    % [1] Theorem 1
    M = ComparisonMatrix(X);
    if IsPSD(M)
        iki = 1;
        return
    elseif k == 2
        iki = 0;
        return
    end

    % [1] (8) (Corollary of Theorem 1)
    if trace(X^2) <= 1/(d-1) % and k > 2
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
        iki = IskIncoh(X, k);
    end

    %%  ISKINCOH    Determines whether or not a quantum state is k-incoherent
    %   This function has two required arguments:
    %     RHO: a mixed quantum state
    %     K: a positive integer
    %   
    %   IKINC = IskCoherent(RHO,K) returns 1 if RHO is K-incoherent and return 0
    %   otherwise. This is checked via semidefinite programming.

    %   requires: CVX (http://cvxr.com/cvx/)
    %   author: Nathaniel Johnston (nathaniel@njohnston.ca), based on code and
    %           conversations with Bartosz Regula
    %   last updated: May 14, 2018

    function ikinc = IskIncoh(rho,k)
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
end