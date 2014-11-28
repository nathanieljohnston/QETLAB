%%  UPBSEPDISTINGUISHABLE    Determines whether or not a UPB is distinguishable by separable measurements
%
%   DIST = UPBSepDistinguishable(U,V,W,...) returns 1 or 0, indicating that
%   the UPB specified by the local vectors U, V, W, ... is or is not
%   distinguishable by separable measurements. U, V, W, ... should be
%   matrices, each with the same number of columns, whose columns are the
%   local vectors of a UPB (i.e., each matrix corresponds to one party).
%
%   See [1] for a description of how this computation is carried out.
%
%   URL: http://www.qetlab.com/UPBSepDistinguishable
%
%   References:
%   [1] S. Bandyopadhyay, A. Cosentino, N. Johnston, V. Russo, J. Watrous,
%       and N. Yu. Limitations on separable measurements by convex
%       optimization. E-print: arXiv:1408.6981 [quant-ph], 2014.

%   requires: cvx (http://cvxr.com/cvx/), opt_args.m, vec_partitions.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   last updated: October 30, 2014

function dist = UPBSepDistinguishable(varargin)

    % First, compute all of the "replacement vectors" for this UPB.
    P = upb_replacement_vectors(varargin{:});

    % Convert the replacement vectors that were found above into the projection
    % onto those vectors.
    ct = 1;
    for j = 1:length(P)
        for k = 1:size(P{j}{1},2)
            Q{ct} = P{j}{1}(:,k);
            for l = 2:length(P{j})
                Q{ct} = kron(Q{ct},P{j}{l}(:,k));
            end
            Q{ct} = Q{ct}*Q{ct}';
            ct = ct + 1;
        end
    end
    dim = length(Q{1});

    % Now solve the LP based on these projections that determines whether or
    % not this UPB is distinguishable by separable measurements.
    cvx_begin sdp quiet
        cvx_precision high;
        variable lam(ct)
        
        P_sum = 0;
        for j = 1:ct-1
            P_sum = P_sum + lam(j)*Q{j};
        end
        
        minimize norm(P_sum - eye(dim),'fro')
        subject to
            lam >= 0;
    cvx_end
    dist = real(cvx_optval);
    
    % Convert the output into a form that the user expects (0 = no, 1 = yes).
    if(dist <= sqrt(eps))
        dist = 1;
    else
        dist = 0;
    end
end

% This function computes all of the "replacement vectors" of a UPB. That
% is, it removes one vector from the UPB and finds all product vectors
% orthogonal to the remaining UPB vectors.
function P = upb_replacement_vectors(varargin)

    p = nargin; % number of parties in the UPB
    s = size(varargin{1},2); % number of states in the UPB
    for j = p:-1:1 % get local dimension of each party
        dim(j) = size(varargin{j},1);
    end

    for j = s:-1:1
        P{j} = cell(p,1);
    end

    % Loop through each of the s states in the UPB, removing them one at a
    % time.
    for j = s:-1:1
        cind = setdiff(1:s,j); % the indices of the states that have not been removed

        pt = vec_partitions(cind,p,dim-1); % all possible partitions of the states among the different parties
        num_pt = length(pt); % number of partitions to loop through

        for k = 1:num_pt % loop through the partitions
            num_orth = 1;
            for l = p:-1:1 % loop through the parties
                V = varargin{l}(:,pt{k}{l}); % these are the local states that we are choosing on the l-th party
                orth_state{l} = null(V.');
                num_orth = num_orth*size(orth_state{l},2);
            end
            
            if(num_orth > 1)
                error('upb_replacement_vectors:NotUPB','The set of local vectors provided do not form a UPB.');
            
            % If the found state really is orthogonal to all others, add it to
            % the output.
            elseif(num_orth == 1)
                % First, make sure that this state isn't already in P{j}.
                if(size(P{j}{1},2) > 0)
                    trow = ones(1,size(P{j}{1},2));
                    for l = p:-1:1
                        trow = trow.*(orth_state{l}'*P{j}{l});
                    end
                    new_col = (max(abs(trow)) < 1 - eps*size(P{j}{1},2)^p);
                else
                    new_col = 1;
                end
                
                % OK, this state really is new.
                if(new_col)
                    for l = p:-1:1
                        P{j}{l} = [P{j}{l},orth_state{l}];
                    end
                end
            end
        end
    end
end