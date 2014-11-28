%%  ISTOTALLYNONSINGULAR    Determines whether or not a matrix is totally nonsingular
%   This function has one required argument:
%     X: a matrix
%
%   ITN = IsTotallyNonsingular(X) is either 1 or 0, indicating that X is or
%   is not totally nonsingular (i.e., all of its square submatrices are
%   nonsingular, within reasonable numerical error).
%
%   This function has two optional input arguments:
%     SUB_SIZES (default 1:min(size(X)))
%     TOL (default max(size(X))*eps(norm(X,'fro')))
%
%   [ITN,WIT] = IsTotallyNonsingular(X,SUB_SIZES,TOL) determines whether or
%   not every r-by-r submatrix of X is nonsingular, where r ranges over all
%   values in the vector SUB_SIZES, and nonsingularity is determined within
%   a tolerance of TOL. If ITN = 0 then WIT is a matrix with two rows and r
%   columns that specifies an r-by-r submatrix of X that is singular.
%
%   URL: http://www.qetlab.com/IsTotallyNonsingular

%   requires: opt_args.m, sporth.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: December 13, 2012

function [itn,wit] = IsTotallyNonsingular(X,varargin)

sX = size(X);
wit = 0;
    
% set optional argument defaults: sub_sizes=1:min(size(X)), tol=max(size(X))*eps(norm(X,'fro'))
[sub_sizes,tol] = opt_args({ 1:min(sX), max(sX)*eps(norm(X,'fro')) },varargin{:});

for j = 1:length(sub_sizes)
    % 1x1 rank can be computed quickly, so do it separately
    if(sub_sizes(j) == 1)
        [r,c] = find(abs(X) < tol,1);
        if(min(size(r)) > 0)
            itn = 0;
            wit = [r;c];
            return;
        end
        
    % larger ranks are slower; just loop on through!
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
                % Using det to test for singularity of a matrix is much
                % faster, but much more error-prone than other functions
                % like rank or cond. Thus we use det in general for speed
                % reasons, and then double-check any results we are unsure
                % of via rank (this gets us the best of both worlds).
                if(abs(det(X(sub_ind_r(kr,:),sub_ind_c(kc,:)))) < tol)
                    [~,rnk] = sporth(X(sub_ind_r(kr,:),sub_ind_c(kc,:)),tol);
                    if(rnk < sub_sizes(j))
                        itn = 0;
                        wit = [sub_ind_r(kr,:);sub_ind_c(kc,:)];
                        return;
                    end
                end
            end
        end
    end
end

itn = 1;