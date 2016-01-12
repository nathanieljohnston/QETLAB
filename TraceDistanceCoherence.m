%%  TraceDistanceCoherence    Computes the trace distance of coherence of a quantum state
%   This function has one required argument:
%     RHO: a pure state vector or a density matrix
%
%   [TDC,D] = TraceDistanceCoherence(RHO) is the trace distance of
%   coherence of the quantum state (density matrix) RHO (i.e., it is the
%   minimal distance in the trace norm between RHO and a diagonal density
%   matrix). The optional output argument D is a vector containing the
%   entries of the diagonal state that is closest to RHO
%   (in other words, TraceNorm(RHO - diag(D)) == TDC).
%
%   URL: http://www.qetlab.com/TraceDistanceCoherence

%   requires: cvx (http://cvxr.com/cvx/)
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: January 12, 2016

function [TDC,D] = TraceDistanceCoherence(rho)

    % Compute the size of RHO. If it's a pure state (vector), we can
    % compute it quickly. If not, we have to use a semidefinite program.
    [m,n] = size(rho);
    ispure = (min(m,n) == 1);

    % If the state is single-qubit, we can compute it faster.
    if(max(m,n) == 1)
        TDC = 0;
        D = 1;
        return
    elseif(m == 2 && n == 2)
        TDC = 2*abs(rho(1,2));
        D = diag(diag(rho));
        return
    else
        if(m == n && rank(rho) == 1)
            [~,~,v] = svd(rho); % v(:,1) is the pure state form of rho
            rho = v(:,1);
            ispure = 1;
        end
        
        % Fore pure states, we can compute it much faster than via SDP.
        % The algorithm is somewhat complicated though, so we contain it in
        % the auxiliary function TraceCoherencePure.
        if(ispure)
            [TDC,D] = TraceCoherencePure(rho);
            return
        end
    end

    % In general, trace distance of coherence is computed by semidefinite
    % programming.
    cvx_begin sdp quiet
        cvx_precision best
        variable d(n)

        minimize norm_nuc(rho - diag(d))

        subject to
            d >= 0;
            sum(d) == 1;
    cvx_end

    D = d/sum(d); % divide by sum(d) for numerical reasons
    TDC = real(cvx_optval);
end

% The following function computes the trace distance coherence of a pure
% state using a method that is much faster than the naive SDP.
function [val,D] = TraceCoherencePure(x)
    x = abs(x(:));
    [xs,ind] = sort(x,'descend'); % sort the entries of x in decending order
    
    n = length(xs);
    
    % Take care of a degenerate case.
    if(n == 1)
        val = 0;
        D = 1;
        return
    end
    
    % Now loop through to try to find k. Use binary search to make this as
    % fast as possible.
    klow = 1;
    khigh = n;
    while 1==1 % loop until we explicitly break out
        k = ceil((khigh+klow)/2);
        s = sum(xs(1:k));
        m = sum(xs((k+1):n).^2);
        r = (s^2-1-k*m);
        q = (r + sqrt(r^2 + 4*k*m*s^2))/(2*k*s);
        
        if(q < xs(k))
            % We have found a potential k!
            % Now compute D.
            D = zeros(n,1);
            for j = 1:k
                D(j) = (xs(j) - q)/(s - k*q);
            end

            D = D(perm_inv(ind));

            % We could compute val = sum(svd(x*x' - D)), but the following
            % method is much quicker.
            val = 2*q/(s - k*q); 

            % We want to increase klow (i.e., look at the top half of
            % our range).
            klow = k;
        end

        if(k == khigh)
            % We're done!
            return
        end
        
        if(q >= xs(k))
            % We want to decrease khigh (i.e., look at the bottom half of
            % our range).
            khigh = k;
        end
    end
end