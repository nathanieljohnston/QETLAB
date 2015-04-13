%%  BELLINEQUALITYMAXQUBITS    Approximates the optimal value of a Bell inequality in qubit (i.e., 2-dimensional quantum) settings
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
%   BMAX = BellInequalityMaxQubits(JOINT_COE,A_COE,B_COE,A_VAL,B_VAL) is an
%   upper bound on the maximum value that the specified Bell inequality can
%   take on in 2-dimensional quantum settings (i.e., when the two parties
%   have access to qubits). This bound is computed using the method of
%   arXiv:1308.3410
%
%   This function has now optional input arguments.
%
%   URL: http://www.qetlab.com/BellInequalityMaxQubits

%   requires: CVX (http://cvxr.com/cvx/), PartialTranspose.m,
%             PermutationOperator.m, Swap.m
%
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: April 13, 2015

function [bmax,rho] = BellInequalityMaxQubits(joint_coe,a_coe,b_coe,a_val,b_val)

    % Get some basic values and make sure that the input vectors are column
    % vectors.
    [ma,mb] = size(joint_coe);
    oa = length(a_val);
    ob = length(b_val);
    a_val = a_val(:); b_val = b_val(:);
    a_coe = a_coe(:); b_coe = b_coe(:);
    
    m = ma;
    
    tot_dim = 2^(2*m+2);
    obj_mat = sparse(tot_dim,tot_dim);

    for a = 0:1
        for b = 0:1
            for x = 1:m
                for y = 1:m
                    b_coeff = joint_coe(x,y)*a_val(a+1)*b_val(b+1);
                    if(y==1)
                        b_coeff = b_coeff + a_coe(x)*a_val(a+1);
                    end
                    if(x==1)
                        b_coeff = b_coeff + b_coe(y)*b_val(b+1);
                    end
                    obj_mat = obj_mat + b_coeff*kron(MN_matrix(m,a,x),MN_matrix(m,b,y));
                end
            end
        end
    end
    obj_mat = (obj_mat+obj_mat')/2; % avoids some numerical problems in CVX
    aux_mat = [1,0,0,0; 0,0,0,0; 0,0,0,0; 0,0,0,0];
    
    % Now construct the SDP that does the separability (really PPT)
    % optimization
    cvx_begin sdp quiet
        variable W(2^(2*m),2^(2*m)) hermitian;

        maximize trace(Swap(Swap(kron(Swap(W,[2,m+1],2*ones(1,2*m)),aux_mat),[m+1,2*m+1],2*ones(1,2*m+2)),[m+2,2*m+1],2*ones(1,2*m+2)) * obj_mat)

        subject to
            trace(W) == 1;
            W >= 0;
            
            % Now construct all of the different possible PPT constraints
            for sz = 1:(m-1)
                ppt_partitions = nchoosek(1:(2*m-1),sz);
                for j = 1:size(ppt_partitions,1)
                    PartialTranspose(W, ppt_partitions(j,:), [4,2*ones(1,2*(m-1))]) >= 0;
                end
            end
    cvx_end
    rho = W;
        
    bmax = cvx_optval;
        
    % Deal with error messages.
    if(strcmpi(cvx_status,'Inaccurate/Solved'))
        warning('BellInequalityMaxQubits:NumericalProblems','Minor numerical problems encountered by CVX. Consider adjusting the tolerance level TOL and re-running the script.');
    elseif(strcmpi(cvx_status,'Inaccurate/Infeasible'))
        warning('BellInequalityMaxQubits:NumericalProblems','Minor numerical problems encountered by CVX. Consider adjusting the tolerance level TOL and re-running the script.');
    elseif(strcmpi(cvx_status,'Unbounded') || strcmpi(cvx_status,'Inaccurate/Unbounded') || strcmpi(cvx_status,'Failed'))
        error('BellInequalityMaxQubits:NumericalProblems',strcat('Numerical problems encountered (CVX status: ',cvx_status,'). Please try adjusting the tolerance level TOL.'));
    end
end

% This function computes the matrices M_a^x and N_b^y. Input arguments:
% M: the number of measurement settings for Alice and Bob
function MN = MN_matrix(m,a,x)
    perm = 1:(m+1);
    perm(1) = x+1;
    perm(x+1) = 1;
    
    MN = a*speye(2^(m+1)) + ((-1)^a)*PermutationOperator(2*ones(1,m+1), perm, 0, 1);
end