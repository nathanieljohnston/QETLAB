%%  CLIQUENUMBER    Bounds the clique number (i.e., maximum size of a clique) of a graph
%   This function has one required input argument:
%     A: the adjacency matrix of a graph
%
%   UB = CliqueNumber(A) is an upper bound on the clique number of the
%        graph with adjacency matrix A. This upper bound is computed via
%        the sum-of-squares hierarchy and thus makes use of semidefinite
%        programming.
%
%   This function has two optional input arguments:
%     K (optional, default 0): a non-negative integer that indicates the
%        level of the hierarchy used to bound the optimal value
%     MODE (optional, default 'sos'): either 'sos' or 'nosdp', indicating
%        whether the upper bound should be computed via the sum-of-squares
%        hierarchy or a no-SDP alternate hierarchy. The 'nosdp' hierarchy
%        is less accurate for a given level of the hierarchy, but is much
%        faster and less memory-intensive to run.
%
%   [UB,LB] = CliqueNumber(A,K,MODE) gives an upper bound (UB) and lower
%             bound (LB) on the clique number of the graph with adjacency
%             matrix A, computed via the K-th level of the hierarchy
%             specified by MODE. Some other easily-computed bounds are
%             incorporated into these bounds too.
%
%   URL: http://www.qetlab.com/CliqueNumber

%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: August 4, 2023

function [ub,lb] = CliqueNumber(A,varargin)
    % set optional argument defaults: K=0, MODE='sos'
    [k,mode] = opt_args({ 0, 'sos' },varargin{:});
    do_sos = strcmpi(mode,'sos');% true means SOS, false means SDP-free alternative

    n = length(A);% number of vertices
    d = 2;% max clique is modelled by a degree-4 homogeneous polynomial, we use d = 4/2 = 2 is half the degree
    p = CopositivePolynomial(A);% the degree-4 polynomial associated with A

    if(nargout > 1)% only compute a lower bound if it was requested, since it takes extra computational time
        if(do_sos)
            [ub,lb] = PolynomialSOS(p,n,d,k,'max');
        else
            [ub,lb] = PolynomialOptimize(p,n,d,k,'max');
        end
        lb = ceil(1/(1-lb) - 0.0000000001);% the 0.0000000001 is there for numerical reasons: we do not want to round to an incorrect max clique bound as a result of numerical imprecision
    else
        if(do_sos)
            ub = PolynomialSOS(p,n,d,k,'max');
        else
            ub = PolynomialOptimize(p,n,d,k,'max');
        end
    end

    if(ub >= 1)
        ub = Inf;
    else
        ub = floor(1/(1-ub) + 0.0000000001);% the 0.0000000001 is there for numerical reasons: we do not want to round to an incorrect max clique bound as a result of numerical imprecision
    end

    % Simple upper bound based on number of edges in the graph
    m = round(sum(sum(A))/2);% number of edges
    ub = min(ub,floor((sqrt(8*m+1)+1)/2 + 0.0000000001));% this formula comes from inverting triangular numbers; the 0.0000000001 is there for numerical reasons

    % Simple lower bounds based on number of edges and other easily-computable properties of the graph
    if(nargout > 1)
        lb = max(lb,ceil(2*m/(2*m - max(abs(eig(A)))^2) - 0.0000000001));% Nikiforov bound
    end
end