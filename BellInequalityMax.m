%%  BELLINEQUALITYMAX    Computes the maximum value of a Bell inequality
%   This function has three required input arguments:
%     COEFFICIENTS: a matrix or tensor that specifies the coefficients of the
%					Bell inequality in either full probability, full correlator,
%					or Collins-Gisin notation
%     DESC: a vector [oa ob ma mb] that specifies the number of outputs and inputs
%			of Alice and Bob
%	  NOTATION: a character vector that specifies the notation the coefficients are
%				written in. the possible values are 'fp', 'fc', and 'cg', for
%				full probability, full correlator, and Collins-Gisin
%
%   BMAX = BellInequalityMax(COEFFICIENTS,DESC,NOTATION) is the
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
%   BMAX = BellInequalityMax(COEFFICIENTS,DESC,NOTATION,MTYPE,K) is
%   the maximum value that the specified Bell inequality can take on in the
%   setting (classical, quantum, or no-signalling) specified by MTYPE.
%
%   URL: http://www.qetlab.com/BellInequalityMax

%   requires: CVX (http://cvxr.com/cvx/)
%   author: Nathaniel Johnston (nathaniel@njohnston.ca) and Mateus Araújo
%   package: QETLAB
%   last updated: September 16, 2022

function bmax = BellInequalityMax(coefficients,desc,notation,varargin)

    % set optional argument defaults: MTYPE='classical', K=1
    [mtype,k] = opt_args({ 'classical', 1 },varargin{:});

    % Get some basic values
    oa = desc(1);
    ob = desc(2);
    ma = desc(3);
    mb = desc(4);
    
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
            
            if (strcmpi(notation,'fp'))
            	M = coefficients;
            elseif (strcmpi(notation,'cg'))
            	M = cg2fp(coefficients,desc);
            elseif (strcmpi(notation,'fc'))
            	M = fc2fp(coefficients);
            end
            
            bmax = M(:)'*p(:);
            
            maximize bmax;

            subject to
                NPAHierarchy(p,k) == 1;
        cvx_end
        
        % Deal with error messages.
        if(strcmpi(cvx_status,'Inaccurate/Solved'))
            warning('BellInequalityMax:NumericalProblems','Minor numerical problems encountered by CVX. Consider adjusting the tolerance level TOL and re-running the script.');
        elseif(strcmpi(cvx_status,'Inaccurate/Infeasible'))
            warning('BellInequalityMax:NumericalProblems','Minor numerical problems encountered by CVX. Consider adjusting the tolerance level TOL and re-running the script.');
        elseif(strcmpi(cvx_status,'Unbounded') || strcmpi(cvx_status,'Inaccurate/Unbounded') || strcmpi(cvx_status,'Failed'))
            error('BellInequalityMax:NumericalProblems',strcat('Numerical problems encountered (CVX status: ',cvx_status,'). Please try adjusting the tolerance level TOL.'));
        end
    elseif(strcmpi(mtype,'classical'))
        % Compute the classical maximum value by going through all of Bob's strategies and picking
        % Alice's corresponding optimal strategy

		if oa == 2 && ob == 2 %if we have only 2 outcomes convert to full correlation notation to use the faster algorithm        
            if (strcmpi(notation,'fc'))
            	M = coefficients;
            elseif (strcmpi(notation,'fp'))
            	M = fp2fc(coefficients);
            elseif (strcmpi(notation,'cg'))
            	M = cg2fc(coefficients);
            end        
        
	        if (ma < mb) %if Alice has fewer inputs than Bob we swap them
		    	M = M';
				ma = desc(4);
				mb = desc(3);
			end
	
			if (mb <= 21) %for few inputs it's faster not to parallelize
		   		parallel_threads = 0;
			else
		   		parallel_threads = Inf;
			end

			Bmarginal = M(1,2:end);
			Amarginal = M(2:end,1);
			correlations = M(2:end,2:end);
			bmax = -Inf;

			parfor (b=0:2^mb-1,parallel_threads)
			    b_vec = 1-2*integer_digits(b,2,mb);
				temp_bmax = Bmarginal*b_vec + sum(abs(Amarginal + correlations*b_vec));
			    bmax = max(bmax,temp_bmax);
			end
			bmax = bmax + M(1,1);

		else
			if (strcmpi(notation,'fp')) %if we have more than 2 output we convert to full probability notation
				M = coefficients;
			elseif (strcmpi(notation,'cg'))
				M = CG2FP(coefficients,desc);
			end
			
			if (oa^ma < ob^mb) % we choose Bob as the party with the fewest strategies
			    M = permute(M,[2,1,4,3]);
			    [oa ob ma mb] = size(M);
			end

			M = permute(M,[1 3 2 4]);
			M = reshape(M,oa*ma,ob*mb);

			offset = 1+ob*[0:mb-1]';

			if (ob^mb <= 10^6) %for few strategies it's faster not to parallelize
			   parallel_threads = 0;
			else
			   parallel_threads = Inf;
			end

			bmax = -Inf;
			parfor (b=0:ob^mb-1,parallel_threads)
		        ind = integer_digits(b,ob,mb) + offset;
		        Ma = sum(M(:,ind),2);
		        temp_bmax = sum(max(reshape(Ma,oa,ma)));
			    bmax = max(bmax,temp_bmax);
			end

		end
	
    else
        error('BellInequalityMax:InvalidMTYPE','MTYPE must be one of ''classical'', ''quantum'', or ''nosignal''.');
    end
end

% converts number "number" to base "base" with digits "digits"
function dits = integer_digits(number,base,digits)
    dits = zeros(digits,1);
    for i=1:digits
        dits(digits+1-i) = mod(number,base);
        number = floor(number/base);
    end
end
