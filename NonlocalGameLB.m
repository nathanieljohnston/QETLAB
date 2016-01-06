%%  NONLOCALGAMELB    Computes a lower bound on the quantum value of a non-local game
%   This function has three required input arguments:
%     D: the local dimension (e.g., D = 2 corresponds to Alice and Bob each
%        having access to a qubit)
%     P: a matrix whose (x,y)-entry is the probability that the referee
%        asks Alice question x and Bob question y
%     V: a 4-D array whose (a,b,x,y)-entry is the value given to Alice and
%        Bob when they provide answers a and b respectively to questions x
%        and y.
%
%   NGLB = NonlocalGameLB(D,P,V) is a lower bound on the maximum value that 
%   the specified non-local game can take on in quantum mechanical settings
%   where Alice and Bob each have access to D-dimensional quantum systems.
%
%   This function works by starting with a randomly-generated POVM for Bob,
%   and then optimizing Alice's POVM and the shared entangled state. Then
%   Alice's POVM and the entangled state are fixed and Bob's POVM is
%   optimized. And so on, back and forth between Alice and Bob until
%   convergence is reached.
%
%   This function has one optional input argument:
%     VERBOSE (default 1): a flag (either 1 or 0) indicating that the
%     function should or should not print out partial progress as it works. 
%
%   URL: http://www.qetlab.com/NonlocalGameLB

%   requires: CVX (http://cvxr.com/cvx/), opt_args.m, opt_disp.m,
%             RandomDensityMatrix.m
%
%   author: Nathaniel Johnston (nathaniel@njohnston.ca) and
%           Vincent Russo (vrusso@uwaterloo.ca), partially based on some
%           code by John Watrous
%   package: QETLAB
%   last updated: June 24, 2015

function nglb = NonlocalGameLB(d,p,V,varargin)

    % set optional argument defaults: VERBOSE=1
    [verbose] = opt_args({ 1 },varargin{:});
    
    % Get some basic values.
    [ma,mb] = size(p);
    oa = size(V,1);
    ob = size(V,2);
    
    % Generate random starting measurements for Bob.    
    B = zeros(d,d,ob,mb);
    for y = 1:mb
        sum_B = zeros(d);
        for b = 1:ob-1
            B(:,:,b,y) = RandomDensityMatrix(d,0,1);
            sum_B = sum_B + B(:,:,b,y);
        end
        lam = norm(sum_B);
        B(:,:,:,y) = B(:,:,:,y)/(lam+0.1); % scale Bob's measurements so that they add up to less than the identity
        B(:,:,ob,y) = eye(d) - sum_B/(lam+0.1); % fill in the final measurement outcome
    end  

    % Now loop until you reach convergence.
    it_diff = 1;
    nglb = -1;
    ct = 0;
    opt_disp('Iter | Value\n-------------\n',verbose);
        
    while it_diff > 10^-6
        % First, optimize over Alice's measurements and the shared entangled state.
        cvx_begin quiet
            variable q(oa,ob,ma,mb);
            variable rho(d,d) hermitian;
            variable A(d,d,oa,ma) hermitian;
            maximize sum(sum(p.*squeeze(sum(sum(V.*q,1),2))))
            subject to
                for a = oa:-1:1
                    for b = ob:-1:1
                        for x = ma:-1:1
                            for y = mb:-1:1
                                q(a,b,x,y) == trace(B(:,:,b,y)'*A(:,:,a,x));
                            end
                        end
                    end
                end
                for a = 1:oa
                    for x = 1:ma
                        A(:,:,a,x) == hermitian_semidefinite(d);
                    end
                end
                A_a_sum = sum(A,3);
                for x = 1:ma
                    A_a_sum(:,:,x) == rho;
                end
                rho == hermitian_semidefinite(d);
                trace(rho) == 1;
        cvx_end
        
        % Next, optimize over Bob's measurements.
        cvx_begin quiet
            variable q(oa,ob,ma,mb);
            variable B(d,d,ob,mb) hermitian;
            maximize sum(sum(p.*squeeze(sum(sum(V.*q,1),2))))
            subject to
                for a = oa:-1:1
                    for b = ob:-1:1
                        for x = ma:-1:1
                            for y = mb:-1:1
                                q(a,b,x,y) == trace(B(:,:,b,y)'*A(:,:,a,x));
                            end
                        end
                    end
                end
                for b = 1:ob
                    for y = 1:mb
                        B(:,:,b,y) == hermitian_semidefinite(d);
                    end
                end
                B_b_sum = sum(B,3);
                for y = 1:mb
                    B_b_sum(:,:,y) == eye(d);
                end
        cvx_end
        
        % Update the best value found so far, so that we can break out if
        % this method has converged.
        it_diff = real(cvx_optval) - nglb;
        nglb = real(cvx_optval);
        ct = ct + 1;
        
        % If requested, display the best bound we have found so far.
        if(ct >= 1000)
            padstr = '';
        elseif(ct >= 100)
            padstr = ' ';
        elseif(ct >= 10)
            padstr = '  ';
        else
            padstr = '   ';
        end
        opt_disp([num2str(ct),padstr,' | ',num2str(nglb,4),'\n'],verbose);
    end
end