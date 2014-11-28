%%  ISTOTALLYPOSITIVE    Determines whether or not a matrix is totally positive
%   This function has one required argument:
%     X: a matrix
%
%   ITP = IsTotallyPositive(X) is either 1 or 0, indicating that X is or
%   is not totally positive (i.e., all of its square submatrices have
%   positive real determinant).
%
%   This function has two optional input arguments:
%     SUB_SIZES (default 1:min(size(X)))
%     TOL (default max(size(X))*eps(norm(X,'fro')))
%
%   [ITP,WIT] = IsTotallyPositive(X,SUB_SIZES,TOL) determines whether or
%   not every r-by-r submatrix of X has positive determinant, where r
%   ranges over all values in the vector SUB_SIZES, and positivity is
%   determined within a tolerance of TOL. If ITP = 0 then WIT is a matrix
%   with two rows and r columns that specifies an r-by-r submatrix of X
%   that does not have positive determinant.
%
%   URL: http://www.qetlab.com/IsTotallyPositive

%   requires: opt_args.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: December 13, 2012

function [itp,wit] = IsTotallyPositive(X,varargin)

sX = size(X);
wit = 0;

% set optional argument defaults: sub_sizes=1:min(size(X)), tol=max(size(X))*eps(norm(X,'fro')) (same as rank)
[sub_sizes,tol] = opt_args({ 1:min(sX), max(sX)*eps(norm(X,'fro')) },varargin{:});

for j = 1:length(sub_sizes)
    % 1x1 determinants can be computed quickly, so do them separately
    if(sub_sizes(j) == 1)
        [r,c] = find(min(real(X),-abs(imag(X))) < -tol,1);
        if(min(size(r)) > 0)
            itp = 0;
            wit = [r;c];
            return;
        end
        
    % larger determinants are slower; just loop on through!
    else
        sub_ind_r = nchoosek(1:sX(1),sub_sizes(j));
        if(sX(1) == sX(2)) % nchoosek is slightly slow, so only call it once if we can get away with it
            sub_ind_c = sub_ind_r;
        else
            sub_ind_c = nchoosek(1:sX(2),sub_sizes(j));
        end
        sub_ind_len_r = size(sub_ind_r,1);
        sub_ind_len_c = size(sub_ind_c,1);
        for kr = 1:sub_ind_len_r
            for kc = 1:sub_ind_len_c
                d = det(X(sub_ind_r(kr,:),sub_ind_c(kc,:)));
                if(d < -tol || abs(imag(d)) > tol)
                    itp = 0;
                    wit = [sub_ind_r(kr,:);sub_ind_c(kc,:)];
                    return;
                end
            end
        end
    end
end

itp = 1;