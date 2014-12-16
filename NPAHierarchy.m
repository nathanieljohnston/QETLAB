%%  NPAHIERARCHY    Determines whether or not a set of probabilities satisfy the conditions of the NPA hierarchy
%   This function has one required input argument:
%     P: a 4D array of probabilities, where P(a,b,x,y) is the probability
%        that Alice and Bob get measurement outcome (a,b) when they use
%        measurement setting (x,y)
%
%   IS_NPA = NPAHierarchy(P) is either 1 or 0, indicating that the
%   probabilities in the array P do or do not satisfy the conditions of the
%   first level of the NPA heierarchy.
%
%   This function has one optional input argument:
%     K (default 1): a non-negative integer
%
%   IS_NPA = NPAHierarchy(P,K) is as above, but the K-th level of the NPA
%   hierarchy is checked. If K = 0 then it is just verified that P is a
%   valid probability array (i.e., each entry is non-negative, each matrix
%   "face" adds up to 1, and Alice's and Bob's marginal probabilities are
%   consistent with each other).
%
%   URL: http://www.qetlab.com/NPAHierarchy

%   requires: CVX (http://cvxr.com/cvx/), opt_args.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: December 9, 2014

function is_npa = NPAHierarchy(p,varargin)

    % set optional argument defaults: K=1
    [k] = opt_args({ 1 },varargin{:});

    % Start by computing the number measurement settings for Alice and Bob (MA
    % and MB) and the number of outcomes for each measurement setting (OA and
    % OB).
    [oa,ob,ma,mb] = size(p);
    o_vec = [oa;ob];
    m_vec = [ma;mb];
    tot_dim = (o_vec'*m_vec)^k; % dimension of the matrix GAMMA used in the NPA SDP
    tol = eps^(3/4); % numerical tolerance used

    % Make sure that P really is a probability array and that its marginal
    % probabilities are consistent with each other (but only if P is not a
    % CVX variable... if it is a CVX variable we will enforce these
    % constraints within the SDP below).
    if(~isa(p,'cvx'))
        % Require that P is a probability matrix (all entries are
        % non-negative and its faces sum to 1).
        if(min(min(min(min(p)))) < -tol^(3/4))
            is_npa = 0;
            return;
        end
        
        for i = 1:ma
            for j = 1:mb
                if(abs(sum(sum(p(:,:,i,j))) - 1) > tol^(3/4))
                    is_npa = 0;
                    return;
                end
            end
        end
        
        % Let's check Bob's marginal probabilities.
        for i = 1:ob
            for j = 1:mb
                marg = sum(squeeze(p(:,i,:,j)),1);
                if(max(abs(marg - marg(1))) > tol^(3/4))
                    is_npa = 0;
                    return;
                end
            end
        end
        
        % Now check Alice's marginal probabilities.
        for i = 1:oa
            for j = 1:ma
                marg = sum(squeeze(p(i,:,j,:)),1);
                if(max(abs(marg - marg(1))) > tol^(3/4))
                    is_npa = 0;
                    return;
                end
            end
        end
    end
    
    % Check the NPA SDP (if K is large enough).
    if(k >= 1 || isa(p,'cvx'))
        i_ind = zeros(2,k);
        j_ind = zeros(2,k);

        cvx_begin quiet
            cvx_precision(tol);
            % We only enforce the actual NPA constraints if K >= 1... if
            % K = 0 we are just verifying marginals are consistent.
            if(k >= 1)
                variable G(tot_dim,tot_dim) symmetric
            end
            subject to
                if(k >= 1)
                    % The following double loop loops over all entries of G and
                    % enforces entry-by-entry the (somewhat complicated) set of NPA
                    % constraints.
                    for i = 1:tot_dim
                        for j = i:tot_dim
                            % First, determine what "type" of product of
                            % measurement operators the given matrix entry
                            % corresponds to (see product_of_orthogonal function
                            % below for details).
                            res = product_of_orthogonal([fliplr(i_ind),j_ind],m_vec);

                            % Entry is 0 if S_i^dagger*S_j = 0.
                            if(isequal(res,0))
                                G(i,j) == 0;

                            % Entry is a single probability from the P array if
                            % S_i^dagger*S_j measures on both Alice and Bob's
                            % sytems.
                            elseif(size(res,2) == 2)
                                res = sortrows(res.').'; % makes sure that Alice's measurement comes first in RES
                                G(i,j) == p(res(2,1)+1,res(2,2)+1,res(1,1)+1,res(1,2)-m_vec(1)+1);

                            % Entry is a sum of probabilities from the P array if
                            % S_i^dagger*S_j measures on just one system.
                            elseif(isequal(size(res),[2,1]))
                                if(res(1) >= m_vec(1)) % measure on Bob's system
                                    G(i,j) == sum(p(:,res(2)+1,1,res(1)-m_vec(1)+1));
                                else % measure on Alice's system
                                    G(i,j) == sum(p(res(2)+1,:,res(1)+1,1));
                                end
                            end

                            j_ind = update_ind(j_ind,k,m_vec,o_vec);
                        end
                        i_ind = update_ind(i_ind,k,m_vec,o_vec);
                        j_ind = i_ind;
                    end
                    G == semidefinite(tot_dim);
                end
                
                % Now enforce that P is a probability array and that its
                % marginals are consistent with each other.
                if(isa(p,'cvx'))
                    % First, P is an array of probabilities, so it must be
                    % non-negative and its faces must sum to 1.
                    p >= 0;
                    for i = 1:ma
                        for j = 1:mb
                            sum(sum(p(:,:,i,j))) == 1;
                        end
                    end
                    
                    % Bob's marginal probabilities must be consistent.
                    for i = ob:-1:1
                        for j = mb:-1:1
                            marg_a{i,j} = sum(squeeze(p(:,i,:,j)),1);
                            marg_a{i,j} == marg_a{i,j}(1)*ones(1,ma);
                        end
                    end

                    % Alice's marginal probabilities must be consistent.
                    for i = oa:-1:1
                        for j = ma:-1:1
                            marg_b{i,j} = sum(squeeze(p(i,:,j,:)),1);
                            marg_b{i,j} == marg_b{i,j}(1)*ones(1,mb);
                        end
                    end
                end
        cvx_end

        is_npa = 1-min(cvx_optval,1); % CVX-safe way to map (0,Inf) to (1,0)
        if(~isa(p,'cvx')) % make the output prettier if it's not a CVX input
            is_npa = round(is_npa);
        end
    else
        is_npa = 1;
    end
end

% This is a function that computes the next index matrix, given an old
% one. Index matrices have 2 rows, and keep track of the measurement
% operators that are being multiplied together, from left to right. The
% first row contains the index of the measurement, and the second row
% contains the index of the result of that measurement. For example,
% the index matrix [2 3 3;0 1 4] refers to the product of 3 measurement
% operators. The first measurement operator is result 0 of measurement
% 2, the second is result 1 of measurement 3, and the third is result 4
% of measurement 3. Note that the entries of this matrix are indexed
% starting at 0 (i.e., there is a measurement 0 and a measurement
% result 0).
%
% Note that we need this function, rather than the simpler update_odometer
% function, because the upper limit in the second row of an index matrix
% depends on the value in the first row.
function new_ind = update_ind(old_ind,k,m_vec,o_vec)
    % Start by increasing the last index by 1.
    new_ind = old_ind;
    new_ind(2,k) = new_ind(2,k)+1;
        
    % Now we work the "odometer": repeatedly set each digit to 0 if it
    % is too high and carry the addition to the left until we hit a
    % digit that *isn't* too high.
    for l = k:-1:1
        % If we've hit the end of the outcomes for this particular
        % measurement, move to the next measurement setting.
        if(new_ind(2,l) >= o_vec(min(floor(new_ind(1,l)/m_vec(1)),1)+1))
            new_ind(2,l) = 0;
            new_ind(1,l) = new_ind(1,l) + 1;
        else
            return; % always return if the odometer doesn't turn over
        end

        % If we've hit the end of the measurement settings within
        % this particular measurement operator, move to the previous
        % measurement operator.
        if(new_ind(1,l) >= sum(m_vec))
            new_ind(1,l) = 0;
            if(l >= 2)
                new_ind(2,l-1) = new_ind(2,l-1) + 1;
            else % if L = 1 too then we are completely done! It doesn't matter what we do from here; just exit.
                return;
            end
        else
            return;
        end
    end
end

% This function determines the nature of the operator specified by the
% index matrix IND. If IND corresponds to something that is not (generally)
% a measurement, then -1 is returned. If it corresponds to the zero
% operator, 0 is returned. If it corresponds to a measurement then the
% index matrix of that measurement is returned.
function res = product_of_orthogonal(ind,m_vec)
    res = -1;
    len = size(ind,2);
    
    % IND is the product of just one measurement operator.
    if(len == 1)
        res = ind;
        return;
        
    % IND is the product of two commuting measurement operators.
    elseif(len == 2 && min(floor(ind(1,1)/m_vec(1)),1) ~= min(floor(ind(1,2)/m_vec(1)),1))
        res = ind;
        return;
    end
    
    % IND is more complicated. Recursively figure out how much it can be
    % simplified.
    for i = 1:len-1
        for j = i+1:len
            % These two measurements are next to each other and are
            % orthogonal!
            if(ind(1,i) == ind(1,j) && ind(2,i) ~= ind(2,j))
                res = 0;
                return; % one is enough; break out of the loop when one is found
                
            % These two measurements are next to each other and are the
            % same! Merge them and then start over.
            elseif(ind(1,i) == ind(1,j) && ind(2,i) == ind(2,j))
                res = product_of_orthogonal(ind(:,[1:j-1,j+1:len]),m_vec);
                return;
                
            % These two measurements act on the same party, but are not
            % orthogonal. Stop increasing J and increase I again.
            elseif(min(floor(ind(1,i)/m_vec(1)),1) == min(floor(ind(1,j)/m_vec(1)),1))
                break;
            end
        end
    end
end