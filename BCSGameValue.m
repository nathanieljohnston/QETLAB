%%  BCSGAMEVALUE    Computes the maximum value of a binary constraint system (BCS) game
%   This function has one required input argument:
%     C: A cell, each of whose elements is a constraint in the BCS. The
%        constraints themselves are specified as 2-by-2-by-2-by-... binary
%        arrays, where the (i,j,k,...)-entry is 1 if and only if setting
%        v1=i, v2=j, v3=k, ... satisfies that constraint.
%
%   BCSVAL = BCSGameValue(C) is the maximum value that the specified 
%   BCS game can take on in classical mechanics. For the maximum quantum
%   or no-signalling value, see the optional arguments described below.
%
%   This function has two optional input arguments:
%     MTYPE (default 'classical'): one of 'classical', 'quantum', or
%       'nosignal', indicating what type of BCS game value should be
%       computed. IMPORTANT NOTE: if MTYPE='quantum' then only an upper
%       bound on the nonlocal game value is computed, not its exact value
%       (see the argument K below).
%     K (default 1): if MYTPE='quantum', then K is a non-negative integer
%       or string indicating what level of the NPA hierarchy to use to
%       bound the nonlocal game (higher values give better bounds, but
%       require more computation time). See the NPAHierarchy function for
%       details.
%
%   BCSVAL = BCSGameValue(C,MTYPE,K) is the maximum value that the
%   specified nonlocal game can take on in the setting (classical, quantum,
%   or no-signalling) specified by MTYPE.
%
%   URL: http://www.qetlab.com/BCSGameValue

%   requires: CVX (http://cvxr.com/cvx/), NonlocalGameValue.m,
%             NPAHierarchy.m, opt_args.m, update_odometer.m
%
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: June 23, 2015

function bcsval = BCSGameValue(C,varargin)

    % set optional argument defaults: MTYPE='classical', K=1
    [mtype,k] = opt_args({ 'classical', 1 },varargin{:});

    % Derive the form of the nonlocal game from the winning conditions of
    % the binary constraint system game.
    num_cons = length(C); % number of constraints
    num_vars = ndims(C{1}); % number of variables
    
    % Determine which variables appear in which constraints
    var_in_cons = ones(num_cons,num_vars);
    for x = 1:num_cons
        for y = num_vars:-1:1
            tindCell{y} = 1:2;
        end
        for y = 1:num_vars
            indCell = tindCell;
            indCell{y} = 1;
            indVec = ones(1,num_vars);
            indVec(y) = 2;

            if(all(C{x} == repmat(C{x}(indCell{:}),indVec))) % winning condition does not depend on the value of this variable
                var_in_cons(x,y) = 0;
            end
        end
    end
    
    % Compute the probability that each given constraint/question pair is
    % asked.
    p = zeros(num_cons,num_vars);
    max_num_vars_in_cons = 0;
    for x = 1:num_cons
        num_vars_in_cons = sum(var_in_cons(x,:));
        max_num_vars_in_cons = max(max_num_vars_in_cons,num_vars_in_cons);
        for y = 1:num_vars
            if(var_in_cons(x,y) == 0) % variable does not appear in constraint
                p(x,y) = 0;
            else
                p(x,y) = 1/(num_cons*num_vars_in_cons);
            end
        end
    end
    
    % Now compute which answer pairs correspond to Alice and Bob winning.
    V = zeros(2^max_num_vars_in_cons,2,num_cons,num_vars); % a = 2^num_vars, b = 2, x = num_cons, y = num_vars
    for a = 1:2^max_num_vars_in_cons
        % Create a cell containing the binary representation of a as its
        % entries, but with 1's and 2's instead of 0's and 1's. Also, pad
        % on the left with 1's so that each of these have the same length
        % as we loop over a.
        a_cell = cellfun(@str2num,num2cell(dec2bin(a-1)));
        a_cell = num2cell(1+padarray(a_cell,[0,max_num_vars_in_cons-length(a_cell)],0,'pre'));
        for b = 1:2
            for x = 1:num_cons
                full_a = expand_a(a_cell,var_in_cons(x,:),num_vars);
                for y = 1:num_vars
                    if(C{x}(full_a{:}) == 1 && b == full_a{y}) % this variable assignment by Alice satisfies constraint x, and this variable assignment by Bob agrees with Alice
                        V(a,b,x,y) = 1;
                    end
                end
            end
        end
    end
    
    % Compute the value of this binary constraint system game, which can be
    % phrased as a general nonlocal game.
    bcsval = NonlocalGameValue(p,V,mtype,k);
end

function full_a = expand_a(small_a,var_list,num_vars)
    ct = 1;
    for j = 1:num_vars
        if(var_list(j) == 1)
            full_a{j} = small_a{ct};
            ct = ct + 1;
        else
            % We don't care about this variable's value; set it to anything.
            full_a{j} = 1;
        end
    end
end