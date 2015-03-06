%%  BELLINEQUALITYMAX    Computes the maximum value of a Bell inequality
%   This function has five required input arguments:
%     JOINT_COE: a matrix whose (i,j)-entry is the coefficient of the term
%                <A_iB_j> in the Bell inequality.
%     A_COE: a vector whose i-th entry is the coefficient of the term <A_i>
%            in the Bell inequality.
%     B_COE: a vector whose i-th entry is the coefficient of the term <B_i>
%            in the Bell inequality.
%     A_VAL: a vector whose i-th entry is the value of the i-th measurement
%            outcome on Alice's system
%     B_VAL: a vector whose i-th entry is the value of the i-th measurement
%            outcome on Bob's system
%
%   BMAX = BellInequalityMax(JOINT_COE,A_COE,B_COE,A_VAL,B_VAL) is the
%   maximum value that the specified Bell inequality can take on in
%   classical mechanics. For the maximum quantum or no-signalling value,
%   see the optional arguments described below.
%
%   This function has two optional input arguments:
%     MTYPE (default 'classical'): one of 'classical', 'quantum', or
%       'nosignal', indicating what type of Bell inequality maximum should
%       be computed. IMPORTANT NOTE: if MTYPE='quantum' then only an upper
%       bound on the Bell inequality is computed, not its exact value (see
%       the argument K below).
%     K (default 1): if MYTPE='quantum', then K is a non-negative integer
%       or string indicating what level of the NPA hierarchy to use to
%       bound the Bell inequality (higher values give better bounds, but
%       require more computation time). See the NPAHierarchy function for
%       details.
%
%   BMAX = BellInequalityMax(JOINT_COE,A_COE,B_COE,A_VAL,B_VAL,MTYPE,K) is
%   the maximum value that the specified Bell inequality can take on in the
%   setting (classical, quantum, or no-signalling) specified by MTYPE.
%
%   URL: http://www.qetlab.com/BellInequalityMax

%   requires: CVX (http://cvxr.com/cvx/), NPAHierarchy.m, opt_args.m,
%             update_odometer.m
%
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: March 6, 2015

function bmax = BellInequalityMax(joint_coe,a_coe,b_coe,a_val,b_val,varargin)

    % set optional argument defaults: MTYPE='classical', K=1
    [mtype,k] = opt_args({ 'classical', 1 },varargin{:});

    % Get some basic values and make sure that the input vectors are column
    % vectors.
    [ma,mb] = size(joint_coe);
    oa = length(a_val);
    ob = length(b_val);
    a_val = a_val(:); b_val = b_val(:);
    a_coe = a_coe(:); b_coe = b_coe(:);
    
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
            variable p(oa,ob,ma,mb);

            % Set up the Bell inequality.
            D = squeeze(sum(sum(repmat(a_val*b_val',[1,1,ma,mb]).*p,1),2));
            Da = sum(repmat(a_val(:),[1,ma]).*squeeze(sum(p(:,:,:,1),2)),1);
            Db = sum(repmat(b_val(:),[1,mb]).*squeeze(sum(p(:,:,1,:),1)),1);

            maximize sum(sum(joint_coe.*D)) + Da*a_coe + Db*b_coe

            subject to
                NPAHierarchy(p,k) == 1;
        cvx_end
        
        bmax = cvx_optval;
        
        % Deal with error messages.
        if(strcmpi(cvx_status,'Inaccurate/Solved'))
            warning('BellInequalityMax:NumericalProblems','Minor numerical problems encountered by CVX. Consider adjusting the tolerance level TOL and re-running the script.');
        elseif(strcmpi(cvx_status,'Inaccurate/Infeasible'))
            warning('BellInequalityMax:NumericalProblems','Minor numerical problems encountered by CVX. Consider adjusting the tolerance level TOL and re-running the script.');
        elseif(strcmpi(cvx_status,'Unbounded') || strcmpi(cvx_status,'Inaccurate/Unbounded') || strcmpi(cvx_status,'Failed'))
            error('BellInequalityMax:NumericalProblems',strcat('Numerical problems encountered (CVX status: ',cvx_status,'). Please try adjusting the tolerance level TOL.'));
        end
    elseif(strcmpi(mtype,'classical'))
        % Compute the classical maximum value just by brute-forcing over
        % all possibilities.
        bmax = -Inf;
        a_ind = zeros(1,ma);
        b_ind = zeros(1,mb);
        
        for i = 1:oa^ma
            for j = 1:ob^mb
                bmax = max(bmax,sum(sum(joint_coe.*(a_val(a_ind+1)*b_val(b_ind+1)'))) + a_val(a_ind+1)'*a_coe + b_val(b_ind+1)'*b_coe);
                b_ind = update_odometer(b_ind, ob*ones(1,mb));
            end
            a_ind = update_odometer(a_ind, oa*ones(1,ma));
        end
    else
        error('BellInequalityMax:InvalidMTYPE','MTYPE must be one of ''classical'', ''quantum'', or ''nosignal''.');
    end
end