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
%     K (default 1): a non-negative integer, or a string indicating which
%       level of the hierarchy to check (see below for details)
%
%   IS_NPA = NPAHierarchy(P,K) is as above, but the K-th level of the NPA
%   hierarchy is checked. If K = 0 then it is just verified that P is a
%   valid probability array (i.e., each entry is non-negative, each matrix
%   "face" adds up to 1, and Alice's and Bob's marginal probabilities are
%   consistent with each other).
%
%   If K is a string, it must be a string of a form like '1+ab+aab', which
%   indicates that the intermediate level of the hierarchy should be used,
%   which uses all products of 1 measurement, all products of one Alice and
%   one Bob measurement, and all products of two Alice and one Bob
%   measurement. Use plus signs to separate the different categories of
%   products, as above. The first character of this string should always be
%   a number, indicating the base level to use.
%
%   URL: http://www.qetlab.com/NPAHierarchy

%   requires: CVX (http://cvxr.com/cvx/), opt_args.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: January 22, 2015

function is_npa = NPAHierarchy(p,varargin)

    % set optional argument defaults: K=1
    [k] = opt_args({ 1 },varargin{:});

    % Parse the input argument K to determine which measurement operators
    % to use. BASE_K is the integer part of the input (i.e., we use all
    % operators that are the product of BASE_K or fewer measurements) and K
    % is the maximum number of products that we ever use (e.g., '1+ab+aab'
    % would result in BASE_K = 1, K = 3).
    if(isnumeric(k))
        base_k = k;
        num_k_compon = 1;
    elseif(isa(k,'char'))
        k_types = textscan(lower(k),'%s','Delimiter','+'); % works with old versions of MATLAB, unlike strsplit
        k_types = k_types{1};
        base_k = str2double(k_types{1});
        
        num_k_compon = length(k_types);
        if(num_k_compon > 1)
            k = base_k;
            for j = 2:num_k_compon
                k_types{j} = strtrim(k_types{j});
                k = max(k,length(k_types{j}));
            end
        else
            k = base_k;
        end
    else
        error('NPAHierarchy:InvalidK','K must be a positive integer or a string.');
    end
    
    % Start by computing the number measurement settings for Alice and Bob (MA
    % and MB) and the number of outcomes for each measurement setting (OA and
    % OB).
    [oa,ob,ma,mb] = size(p);
    o_vec = [oa;ob];
    m_vec = [ma;mb];
    tot_dim = 1 + ((o_vec-1)'*m_vec)^k; % upper bound on the dimension of the (compact version of the) matrix GAMMA used in the NPA SDP
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
        i_ind = [zeros(1,k);-ones(1,k)];
        j_ind = [zeros(1,k);-ones(1,k)];

        if(k >= 1)
            % Start by generating all of the product of measurements that
            % you need.
            ind_catalog = cell(0);
            for j = 1:tot_dim
                [res,res_type] = product_of_orthogonal(j_ind,m_vec);
                res_fnd = find_in_cell(res,ind_catalog);
                                    
                % Make sure that this measurement is (1) new, and (2) valid
                % given the user input.
                if(res_fnd == 0 && res_type ~= 0)
                    is_valid_res = (size(res,2) <= base_k);
                    if(~is_valid_res && num_k_compon >= 2)
                        num_a_res = sum(res(1,:) < m_vec(1));
                        num_b_res = size(res,2) - num_a_res;
                        for i = 2:num_k_compon
                            num_a_fnd = length(find(k_types{i}=='a'));
                            num_b_fnd = length(find(k_types{i}=='b'));
                            if(num_a_res <= num_a_fnd && num_b_res <= num_b_fnd)
                                is_valid_res = true;
                                break;
                            end
                        end
                    end
                    
                    if(is_valid_res)
                        ind_catalog{end+1} = res;
                    end
                end
                        
                j_ind = update_ind(j_ind,k,m_vec,o_vec-1);
            end
            real_dim = length(ind_catalog);
        end
        
        cvx_begin quiet
            cvx_precision(tol);
            % We only enforce the actual NPA constraints if K >= 1... if
            % K = 0 we are just verifying marginals are consistent (i.e.,
            % no-signalling).
            if(k >= 1)
                variable G(real_dim,real_dim) symmetric
            end
            subject to
                if(k >= 1)
                    res_catalog = cell(0);
                    res_loc = cell(0);
                    
                    % The following double loop loops over all entries of G and
                    % enforces entry-by-entry the (somewhat complicated) set of NPA
                    % constraints.
                    for i = 1:real_dim
                        for j = i:real_dim
                            % First determine what "type" of product of
                            % measurement operators the given matrix entry
                            % corresponds to (see product_of_orthogonal function
                            % below for details).
                            [res,res_type] = product_of_orthogonal([fliplr(ind_catalog{i}),ind_catalog{j}],m_vec);
                            
                            % Entry is 0 if S_i^dagger*S_j = 0.
                            if(res_type == 0)
                                G(i,j) == 0;
                                
                            % Entry is a single probability from the P array if
                            % S_i^dagger*S_j measures on both Alice and Bob's
                            % sytems.
                            elseif(res_type == 2)
                                G(i,j) == p(res(2,1)+1,res(2,2)+1,res(1,1)+1,res(1,2)-m_vec(1)+1);

                            % Entry is a sum of probabilities from the P array if
                            % S_i^dagger*S_j measures on just one system.
                            elseif(res_type == 1)
                                if(isequal(res,[0;-1]))
                                    G(i,j) == 1; % identity measurement
                                elseif(res(1) >= m_vec(1)) % measure on Bob's system
                                    G(i,j) == sum(p(:,res(2)+1,1,res(1)-m_vec(1)+1));
                                else % measure on Alice's system
                                    G(i,j) == sum(p(res(2)+1,:,res(1)+1,1));
                                end

                            % Entry is a product of non-commuting
                            % measurement operators. We can't specify its
                            % value, but we can specify that it is equal to
                            % other entries that are the *same* product of
                            % measurement operators.
                            else % res_type == -1
                                % Check to see if we have run into this
                                % particular RES before.
                                res_fnd = find_in_cell(res,res_catalog);
                                if(res_fnd == 0)
                                    res_fnd = find_in_cell(product_of_orthogonal(fliplr(res),m_vec),res_catalog);
                                end

                                % No, this RES is new to us.
                                if(res_fnd == 0)
                                    res_catalog{end+1} = res;
                                    res_loc{end+1} = [i,j];

                                % Yes, we have seen this RES before.
                                else
                                    G(i,j) == G(res_loc{res_fnd}(1),res_loc{res_fnd}(2));
                                end
                            end
                        end
                    end
                    G == semidefinite(real_dim);
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
            
            % Deal with error messages.
            if(strcmpi(cvx_status,'Inaccurate/Solved'))
                warning('NPAHierarchy:NumericalProblems','Minor numerical problems encountered by CVX. Consider adjusting the tolerance level TOL and re-running the script.');
            elseif(strcmpi(cvx_status,'Inaccurate/Infeasible'))
                warning('NPAHierarchy:NumericalProblems','Minor numerical problems encountered by CVX. Consider adjusting the tolerance level TOL and re-running the script.');
            elseif(strcmpi(cvx_status,'Unbounded') || strcmpi(cvx_status,'Inaccurate/Unbounded') || strcmpi(cvx_status,'Failed'))
                error('NPAHierarchy:NumericalProblems',strcat('Numerical problems encountered (CVX status: ',cvx_status,'). Please try adjusting the tolerance level TOL.'));
            end
        end
    else
        is_npa = 1;
    end
end


function ind = find_in_cell(val,A)
    ind = 0;
    for i = 1:numel(A)
        if(isequal(A{i},val))
            ind = i;
            return;
        end
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
    % Do we have the identity measurement right now? Go to the first
    % non-identity one.
    if(min(min(old_ind)) == -1)
        new_ind = zeros(2,k);
        return;
    end
    
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
% operator, 0 is returned. If it corresponds to the identity operator,
% [0;0] is returned. If it corresponds to a measurement then the index
% matrix of that measurement is returned.
function [res,res_type] = product_of_orthogonal(ind,m_vec)
    res_type = -1;
    res = ind;
    len = size(ind,2);
    
    % IND is the product of just one measurement operator.
    if(len == 1)
        res_type = 1;
        return;
        
    % IND is the product of two commuting non-identity measurement
    % operators.
    elseif(len == 2 && ind(2,1) >= 0 && ind(2,2) >= 0 && min(floor(ind(1,1)/m_vec(1)),1) ~= min(floor(ind(1,2)/m_vec(1)),1))
        res = sortrows(res.').'; % sort so that Alice's measurement comes first
        res_type = 2;
        return;
    end
    
    % IND is more complicated. Recursively figure out how much it can be
    % simplified.
    for i = 1:len-1
        for j = i+1:len
            % These two measurements are next to each other and are
            % orthogonal!
            if(ind(2,i) >= 0 && ind(2,j) >= 0 && ind(1,i) == ind(1,j) && ind(2,i) ~= ind(2,j))
                res_type = 0;
                return; % one is enough; break out of the loop when one is found

            % These two measurements are next to each other and are the
            % same! Merge them and then start over.
            elseif(ind(1,i) == ind(1,j) && ind(2,i) == ind(2,j))
                [res,res_type] = product_of_orthogonal(ind(:,[1:j-1,j+1:len]),m_vec);
                return;

            % The first of these two measurement operators is the identity.
            % Merge them and then start over.
            elseif(ind(2,i) == -1)
                [res,res_type] = product_of_orthogonal(ind(:,[1:i-1,i+1:len]),m_vec);
                return;

            % The second of these two measurement operators is the
            % identity. Merge them and then start over.
            elseif(ind(2,j) == -1)
                [res,res_type] = product_of_orthogonal(ind(:,[1:j-1,j+1:len]),m_vec);
                return;

            % These two measurements act on the same party, but are not
            % orthogonal. Stop increasing J and increase I again.
            elseif(min(floor(ind(1,i)/m_vec(1)),1) == min(floor(ind(1,j)/m_vec(1)),1))
                break;

            % These two measurements act on different parties and are both
            % non-identity. They commute; move Alice's before Bob's
            elseif(ind(1,i) > ind(1,j))
                [res,res_type] = product_of_orthogonal(ind(:,[1:i-1,j,i+1:j-1,i,j+1:len]),m_vec);
                return;
            end
        end
    end
end