%%  WERNERSTATE    Produces a Werner state
%   This function has two required arguments:
%     DIM: the local dimension
%     ALPHA: the parameters of the Werner state
%
%   RHO = WernerState(DIM,ALPHA) is the Werner state with parameter ALPHA
%   acting on (DIM*DIM)-dimensional space. More specifically, RHO is the
%   density operator defined by (I - ALPHA*S) (normalized to have trace 1),
%   where I is the identity operator and S is the operator that swaps two
%   copies of DIM-dimensional space (see Swap.m and SwapOperator.m for
%   example).
%
%   If ALPHA is a vector with p!-1 entries for some integer p>1 then a
%   multipartite Werner state is returned. This multipartite Werner state
%   is the normalization of I - ALPHA(1)*P(2) - ... - ALPHA(p!-1)*P(p!),
%   where P(i) is the operator that permutes p subsystems according to the
%   i-th permutation when they are written in lexicographical order (for
%   example, the lexicographical ordering when p = 3 is: [1 2 3], [1 3 2],
%   [2 1 3], [2 3 1], [3 1 2], [3 2 1], so P(4) in this case equals
%   PermutationOperator(DIM,[2 3 1]))
%
%   URL: http://www.qetlab.com/WernerState

%   requires: iden.m, opt_args.m, PermutationOperator.m, PermuteSystems.m,
%             Swap.m, SwapOperator.m
%
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: November 12, 2014

function rho = WernerState(dim,alpha)

% compute the Werner state
n_fac = length(alpha)+1; % total number of permutation operators
if(n_fac > 2) % multipartite Werner state
    % compute the number of parties from length(ALPHA)
    n = n_fac;
    for j = 2:n_fac % we won't actually go all the way to n_fac
        n = n/j;
        if(n == j+1)
            break;
        elseif(n < j)
            error('WernerState:InvalidAlpha','The ALPHA vector must contain p!-1 entries for some integer p>1.');
        end
    end
    
    % Done error checking and computing the number of parties -- now
    % compute the Werner state.
    p = sortrows(perms(1:n));
    for j = 2:n_fac
        rho = speye(dim^n) - alpha(j-1)*PermutationOperator(dim,p(j,:),0,1);
    end
    rho = rho/trace(rho);
else % bipartite Werner state
    rho = (speye(dim^2) - alpha*SwapOperator(dim,1))/(dim*(dim-alpha));
end