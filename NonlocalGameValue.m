%%  NONLOCALGAMEVALUE    Computes the maximum value of a two-player nonlocal game
%   This function has two required input arguments:
%     P: a matrix whose (x,y)-entry is the probability that the referee
%        asks Alice question x and Bob question y
%     V: a 4-D array whose (a,b,x,y)-entry is the value given to Alice and
%        Bob when they provide answers a and b respectively to questions x
%        and y.
%
%   NGVAL = NonlocalGameValue(P,V) is the maximum value that the specified 
%   nonlocal game can take on in classical mechanics. For the maximum
%   quantum or no-signalling value, see the optional arguments described
%   below.
%
%   This function has two optional input arguments:
%     MTYPE (default 'classical'): one of 'classical', 'quantum', or
%       'nosignal', indicating what type of nonlocal game value should be
%       computed. IMPORTANT NOTE: if MTYPE='quantum' then only an upper
%       bound on the nonlocal game value is computed, not its exact value
%       (see the argument K below).
%     K (default 1): if MYTPE='quantum', then K is a non-negative integer
%       or string indicating what level of the NPA hierarchy to use to
%       bound the nonlocal game (higher values give better bounds, but
%       require more computation time). See the NPAHierarchy function for
%       details.
%
%   NGVAL = NonlocalGameValue(P,V,MTYPE,K) is the maximum value that the
%   specified nonlocal game can take on in the setting (classical, quantum,
%   or no-signalling) specified by MTYPE.
%
%   URL: http://www.qetlab.com/NonlocalGameValue

%   requires: CVX (http://cvxr.com/cvx/), NPAHierarchy.m, opt_args.m,
%             update_odometer.m
%
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: March 6, 2015

function ngval = NonlocalGameValue(p,V,varargin)

    % set optional argument defaults: MTYPE='classical', K=1
    [mtype,k] = opt_args({ 'classical', 1 },varargin{:});

    % Get some basic values.
    [ma,mb] = size(p);
    oa = size(V,1);
    ob = size(V,2);
    
    % The no-signalling maximum is just implemented by taking the zero-th
    % level of the NPA (quantum) hierarchy.
    if(strcmpi(mtype,'nosignal'))
        mtype = 'quantum';
        k = 0;
    end
    
    % Compute the maximum value of the Bell inequality, depending on which
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
        % Compute the classical maximum value just by brute-forcing over
        % all possibilities.
        ngval = -Inf;
        a_ind = zeros(1,ma);
        b_ind = zeros(1,mb);
        a_ind(ma)=oa-2;
        b_ind(mb)=ob-2;
        p = permute(repmat(p,[1,1,oa,ob]),[3,4,1,2]);
        ab_prob = zeros(oa,ob,ma,mb);
        
        for i = 1:oa^ma
            for j = 1:ob^mb
                for x = 1:ma
                    for y = 1:mb
                        ab_prob(:,:,x,y) = zeros(oa,ob);
                        ab_prob(a_ind(x)+1,b_ind(y)+1,x,y) = 1;
                    end
                end
                tgval = sum(sum(sum(sum(p.*V.*ab_prob))));
                if(tgval >= ngval - 0.000001)
                    ngval = max(ngval, sum(sum(sum(sum(p.*V.*ab_prob)))));
                end
                b_ind = update_odometer(b_ind, ob*ones(1,mb));
            end
            a_ind = update_odometer(a_ind, oa*ones(1,ma));
        end
    else
        error('NonlocalGameValue:InvalidMTYPE','MTYPE must be one of ''classical'', ''quantum'', or ''nosignal''.');
    end
end