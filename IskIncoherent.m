%%  ISKINCOHERENT    Determines whether or not a quantum state is k-incoherent
%   This function has two required arguments:
%     RHO: a mixed quantum state
%     K: a positive integer
%   
%   IKINC = IskCoherent(RHO,K) returns 1 if RHO is K-incoherent and return 0
%   otherwise. This is checked via semidefinite programming.

%   requires: CVX (http://cvxr.com/cvx/)
%   author: Nathaniel Johnston (nathaniel@njohnston.ca), based on code and
%           conversations with Bartosz Regula
%   package: QETLAB
%   last updated: May 14, 2018

function ikinc = IskIncoherent(rho,k)
    n = size(rho,1);

    Pk = nchoosek(1:n,k);
    s = size(Pk,1);

    cvx_begin sdp quiet
    cvx_precision best;
    variable A(k,k,s) hermitian;
    subject to
        P = zeros(n);
        for j=1:s
            proj = zeros(k,n);
            for i = 1:k
                proj(i,Pk(j,i)) = 1;
            end
            P = P + proj'*A(:,:,j)*proj;
            
            A(:,:,j) >= 0;
        end
        rho == P;
    cvx_end
    
    ikinc = 1-min(cvx_optval,1);
end