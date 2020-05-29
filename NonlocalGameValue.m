%%  NONLOCALGAMEVALUE    Computes the maximum value of a two-player non-local game
%   This function has two required input arguments:
%     P: a matrix whose (x,y)-entry is the probability that the referee
%        asks Alice question x and Bob question y
%     V: a 4-D array whose (a,b,x,y)-entry is the value given to Alice and
%        Bob when they provide answers a and b respectively to questions x
%        and y.
%
%   NGVAL = NonlocalGameValue(P,V) is the maximum value that the specified 
%   non-local game can take on in classical mechanics. For the maximum
%   quantum or no-signalling value, see the optional arguments described
%   below.
%
%   This function has three optional input arguments:
%     MTYPE (default 'classical'): one of 'classical', 'quantum', or
%       'nosignal', indicating what type of nonlocal game value should be
%       computed. IMPORTANT NOTE: if MTYPE='quantum' then only an upper
%       bound on the non-local game value is computed, not its exact value
%       (see the argument K below).
%     K (default 1): if MYTPE='quantum', then K is a non-negative integer
%       or string indicating what level of the NPA hierarchy to use to
%       bound the non-local game (higher values give better bounds, but
%       require more computation time). See the NPAHierarchy function for
%       details. If MTYPE is anything else, then K is simply ignored.
%     REPT (default 1): the number of times that the non-local game should
%       be played in parallel (i.e., the number of repetitions of the
%       game).
%
%   NGVAL = NonlocalGameValue(P,V,MTYPE,K,REPT) is the maximum value that
%   the specified non-local game can take on in the setting (classical,
%   quantum or no-signalling) specified by MTYPE.
%
%   URL: http://www.qetlab.com/NonlocalGameValue

%   requires: CVX (http://cvxr.com/cvx/), NPAHierarchy.m, opt_args.m,
%             update_odometer.m
%
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%           Mateus Araújo (wrote classical value computation)
%   package: QETLAB
%   last updated: May 29, 2020

function ngval = NonlocalGameValue(p,V,varargin)

    % set optional argument defaults: MTYPE='classical', K=1, REPT=1
    [mtype,k,rept] = opt_args({ 'classical', 1, 1 },varargin{:});

    % Get some basic values.
    [ma,mb] = size(p);
    oa = size(V,1);
    ob = size(V,2);
    
    % Are we playing the game multiple times? If so, adjust p and V
    % accordingly.
    if(rept > 1)
        % Create the new probability array.
        p = Tensor(p,rept);
        
        % Create the new winning output array (this is more complicated
        % since MATLAB doesn't have any built-in functions for tensoring
        % together 4D arrays).
        V2 = zeros(oa^rept,ob^rept,ma^rept,mb^rept);
        i_ind = zeros(1,rept);
        j_ind = zeros(1,rept);
        for i = 1:ma^rept
            for j = 1:mb^rept
                for l = rept:-1:1
                    to_tensor{l} = V(:,:,i_ind(l)+1,j_ind(l)+1);
                end
                V2(:,:,i,j) = Tensor(to_tensor);

                j_ind = update_odometer(j_ind,mb*ones(1,rept));
            end
            i_ind = update_odometer(i_ind,ma*ones(1,rept));
        end
        V = V2;

        % Update some computed values.
        ma = ma^rept;
        mb = mb^rept;
        oa = oa^rept;
        ob = ob^rept;
    end
    
    % The no-signalling maximum is just implemented by taking the zero-th
    % level of the NPA (quantum) hierarchy.
    if(strcmpi(mtype,'nosignal'))
        mtype = 'quantum';
        k = 0;
    end
    
    % Compute the maximum value of the non-local game, depending on which
    % type of maximum was requested.
    if(strcmpi(mtype,'quantum'))
        cvx_begin quiet
            variable q(oa,ob,ma,mb);
            maximize sum(sum(p.*squeeze(sum(sum(V.*q,1),2))))
            subject to
                NPAHierarchy(q,k) == 1;
        cvx_end
        
        ngval = cvx_optval;
        
        % Deal with error messages.
        if(strcmpi(cvx_status,'Inaccurate/Solved'))
            warning('NonlocalGameValue:NumericalProblems','Minor numerical problems encountered by CVX. Consider adjusting the tolerance level TOL and re-running the script.');
        elseif(strcmpi(cvx_status,'Inaccurate/Infeasible'))
            warning('NonlocalGameValue:NumericalProblems','Minor numerical problems encountered by CVX. Consider adjusting the tolerance level TOL and re-running the script.');
        elseif(strcmpi(cvx_status,'Unbounded') || strcmpi(cvx_status,'Inaccurate/Unbounded') || strcmpi(cvx_status,'Failed'))
            error('NonlocalGameValue:NumericalProblems',strcat('Numerical problems encountered (CVX status: ',cvx_status,'). Please try adjusting the tolerance level TOL.'));
        end
    elseif(strcmpi(mtype,'classical'))
        % Compute the classical maximum value just by looping over Bob's strategies
        % and taking Alice's optimal strategy for each.
        for x=1:ma
            for y=1:mb
                V(:,:,x,y) = p(x,y)*V(:,:,x,y);
            end
        end
        ngval = -Inf;
        b_ind = zeros(1,mb);
        if (oa^ma < ob^mb) % we choose Bob as the party with the fewest strategies
            V = permute(V,[3,4,1,2]);
            [oa,ob,ma,mb] = size(V);
        end
        V = permute(V,[1 3 2 4]); % this is necessary to avoid having to call squeeze later

        if (ob^mb <= 10^6)  % for small problems we don't parallelize
           parallel_threads = 0;
        else
           parallel_threads = Inf;
        end

        parfor (i = 0:ob^mb-1, parallel_threads)
            b_ind = integer_digits(i,ob,mb);
            Valice = zeros(oa,ma);
            for y=1:mb
                Valice = Valice + V(:,:,b_ind(y)+1,y);
            end
            tgval = sum(max(Valice));
            ngval = max(ngval, tgval);
        end
    else
        error('NonlocalGameValue:InvalidMTYPE','MTYPE must be one of ''classical'', ''quantum'', or ''nosignal''.');
    end
end

% converts number "number" to base "base" with digits "digits"
function dits = integer_digits(number,base,digits)
    dits = zeros(1,digits);
    for i=1:digits
        dits(digits+1-i) = mod(number,base);
        number = floor(number/base);
    end
end