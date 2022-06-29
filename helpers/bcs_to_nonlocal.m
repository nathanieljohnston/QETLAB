%%  BCS_TO_NONLOCAL    Converts a constraint description of a binary constraint system (BCS) game to a general non-local game description
%   This function has one required input argument:
%     C: A cell, each of whose elements is a constraint in the BCS. The
%        constraints themselves are specified as 2-by-2-by-2-by-... binary
%        arrays, where the (i,j,k,...)-entry is 1 if and only if setting
%        v1=i, v2=j, v3=k, ... satisfies that constraint.
%
%   [P,V] = BCS_TO_NONLOCAL(C) is the general non-local game description of
%   the BCS game described by the constraints C. That is, P and V can be
%   used as the input to functions like NonlocalGameValue.
%
%   URL: http://www.qetlab.com/bcs_to_nonlocal

%   requires: update_odometer.m
%
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: June 24, 2015

function [p,V] = bcs_to_nonlocal(C)

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
        a_cell = num2cell(1+pad_array(a_cell,max_num_vars_in_cons-length(a_cell), 0));
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