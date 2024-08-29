%%  ISABSKINCOH    Determines whether or not a quantum state is absolutely k-incoherent
%   This function has two required arguments:
%     X: a density matrix
%     k: a positive integer, the absolute coherence level to check for
%
%   IAKI = ISABSKINCOH(X, K) returns 1 if X is absolutely k-incoherent,
%   0 if X is not absolutely k-incoherent, and -1 if the function cannot
%   determine if X is absolutely k-incoherent or not.
%   
%   URL: http://www.qetlab.com/IsAbskIncoh
%   
%   References:
%   [1] N. Johnston, S. Moein, R. Pereira, and S. Plosker. Absolutely 
%       k-Incoherent Quantum States and Spectral Inequalities for Factor Width of 
%       a Matrix. Physical Review A, 106:052417, 2022.

%   requires: CVX (cvxr.com/cvx/)
%   author: Benjamin Talbot, with code from Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: August 26, 2024

function iaki = IsAbskIncoh(X, k)
    iaki = -1;
    tol = 1e-15;
    
    % first check for valid quantum state
    if ~(IsPSD(X) && (trace(X) == 1))
        iaki = 0;
    end

    % trivial: every density matrix is absolutely n-incoherent
    n = size(X, 1);
    if k == n
        iaki = 1;
        return
    end

    eigvals = eig(X);
    rankX = rank(X);
    lmax = max(eigvals);

    % trivial: only the maximally mixed state is absolutely 1-incoherent
    if k == 1
        if all(abs(eigvals - ones(n, 1)/n) <= tol)
            iaki = 1;
            return
        else
            iaki = 0;
            return
        end
    end

    % [1] Theorem 4
    if rankX <= n - k
        iaki = 0;
        return
    elseif rankX == n - k + 1
        if eqnonzeigvals(eigvals)
            iaki = 1;
            return
        end
    end

    % [1] Theorem 5
    if lmax <= 1/(n-k+1)
        iaki = 1;
        return
    end

    
    if k == 2
        % [1] Theorem 7
        if norm(X,'fro')^2 <= 1/(n-1)
            iaki = 1;
            return
        else
            if n <= 3
                iaki = 0;
                return
            end
        end
    elseif k == n - 1
        % [1] Corollary 1
        if lmax > 1 - 1/n
            iaki = 0;
            return
        else
            % [1] Theorem 8
            iaki = IsAbsn1Incoh(eigvals);
        end
    end

    %%  IsAbsn1Incoh    Determines whether or not a vector of eigenvalues corresponds to an absolutely (n-1)-incoherent matrix
    %   This function has one required argument:
    %     LAM: A vector of eigenvalues of the quantum state being tested.
    %          Assumed to be real, but does not have to be sorted.
    %
    %   Returns either 1 or 0, indicating "yes" or "no".
    %
    %   See paper "Absolutely k-Incoherent Quantum States and Spectral
    %   Inequalities for Factor Width of a Matrix" by Johnston et al for
    %   details.

    %   requires: CVX (cvxr.com/cvx/)
    %   author: Nathaniel Johnston (nathaniel@njohnston.ca)
    %   last updated: April 28, 2022

    function iai = IsAbsn1Incoh(lam)
        n = length(lam);
        lam = sort(real(lam),'descend');

        cvx_begin sdp quiet
            cvx_precision best;
            variable L(n,n) symmetric;
            minimize 1;

            subject to
                L(1,1) == -lam(1) - sum(L(1,2:n)) - sum(L(2:n,1));
                for j = 2:n
                    L(j,j) == lam(j);
                end
                L >= 0;
        cvx_end
        
        if(cvx_optval == 1)
            iai = 1;
        else
            iai = 0;
        end
    end

    % Helper function to check if all non-zero eigenvalues are equal
    function allequal = eqnonzeigvals(lam)
        prevlam = lam(1);
        i = 2;
        while prevlam < tol
            prevlam = lam(i);
            i = i + 1;
        end
        allequal = 1;
        for j = i:length(lam)
            if lam(j) > tol % if lam(j) ~= 0
                if lam(j) - prevlam > tol
                    allequal = 0;
                    break
                else
                    prevlam = lam(j);
                end
            end
        end
    end
end