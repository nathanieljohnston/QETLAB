%%  SUPEROPERATOR_DIMS    Computes the input, output, and environment dimensions of a superoperator
%   This function has one required input argument:
%     PHI: a superoperator, represented as either a Choi matrix or a cell
%          of Kraus operators
%
%   This function has three output argument:
%     DA: the dimensions of the input (Alice's) space
%     DB: the dimensions of the output (Bob's) space
%     DE: the dimension of the environment (Eve's) space
%
%   [DA,DB,DE] = superoperator_dims(PHI) returns the dimensions of the
%   input, output, and environment spaces of PHI, in that order. DA and DB
%   are both 1-by-2 vectors containing the row and column dimensions of
%   their spaces. DE is always a scalar, and it is equal to the number of
%   Kraus operators of PHI (if PHI is provided as a Choi matrix then DE is
%   the *minimal* number of Kraus operators of any representation of PHI).
%
%   This function has two optional input arguments:
%     ALLOW_RECT: a flag (either 1 or 0) indicating that the input and
%                 output spaces of PHI can be non-square (default 1)
%     DIM: a vector or matrix containing the input and output dimensions of
%          PHI
%
%   [DA,DB,DE] = superoperator_dims(PHI,ALLOW_RECT,DIM) is as above. DIM
%   should provided if and only if PHI is a Choi matrix with unequal input
%   and output dimensions (since it is impossible to determine the input
%   and output dimensions from the Choi matrix alone). If ALLOW_RECT == 0
%   and PHI acts on non-square matrix spaces, an error will be produced. If
%   PHI maps M_{r,c} to M_{x,y} then DIM should be the 2-by-2 matrix
%   [r,x;c,y]. If PHI maps M_m to M_n, then DIM can simply be the vector
%   [m,n]. If ALLOW_RECT == 0 then DA and DB will be scalars instead of
%   vectors.
%
%   URL: http://www.qetlab.com/superoperator_dims

%   requires: opt_args.m, sporth.m
%   authors: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: November 27, 2014

function [da,db,de] = superoperator_dims(Phi,varargin)

    % If PHI is a cell, the dimensions are easy to compute.
    if(iscell(Phi))
        % Start by computing row and environment dimensions.
        [db,da] = size(Phi{1,1});
        [de,phi_cols] = size(Phi);

        % We now have the row dimensions, so compute the column dimensions.
        if(phi_cols > 2)
            error('superoperator_dims:InvalidDims','The cell PHI should have at most 2 columns (corresponding to left and right Kraus operators).');
        elseif(phi_cols == 2)
            [db(2),da(2)] = size(Phi{1,2});
        else
            da = [da,da];
            db = [db,db];
        end

        % Set optional argument defaults: allow_rect=1, dim=[da',db']
        [allow_rect,dim] = opt_args({ 1, [da',db'] },varargin{:}); % Did the user provide dimensions? Get them.
        dim = expand_dim(dim);
        
        % Now do some error checking.
        if((da(1) ~= da(2) || db(1) ~= db(2)) && allow_rect ~= 1)
            error('superoperator_dims:InvalidDims','The input and output spaces of PHI must be square.');
        elseif(dim(1,1) ~= da(1) || dim(2,1) ~= da(2) || dim(1,2) ~= db(1) || dim(2,2) ~= db(2))
            error('superoperator_dims:InvalidDims','The dimensions of PHI do not match those provided in the DIM argument.');
        end
        for j = 1:de
            if(~all(size(Phi{j,1}) == [dim(1,2),dim(1,1)]) || (phi_cols == 2 && ~all(size(Phi{j,2}) == [dim(2,2),dim(2,1)])))
                error('superoperator_dims:InvalidDims','The Kraus operators of PHI do not all have the same size.');
            end
        end

    % If Phi is a Choi matrix, the dimensions are a bit more of a pain: we have
    % to guess a bit if the input and output dimensions are different.
    else
        % Try to guess DA and DB.
        [r,c] = size(Phi);
        da = [round(sqrt(r)),round(sqrt(c))];
        db = da;

        % Set optional argument defaults: allow_rect=1, dim=[da',db']
        [allow_rect,dim] = opt_args({ 1, [da',db'] },varargin{:});
        dim = expand_dim(dim);
        
        if(r ~= c && allow_rect ~= 1)
            error('superoperator_dims:InvalidDims','The Choi matrix of PHI must be square.');
        end
        if(dim(1,1)*dim(1,2) ~= r || dim(2,1)*dim(2,2) ~= c)
            error('superoperator_dims:InvalidDims','If the input and output dimensions are unequal and PHI is provided as a Choi matrix, the optional argument DIM must be specified (and its dimensions must agree with PHI).');
        end

        [~,de] = sporth(Phi); % environment dimension is the rank of the Choi matrix
    end

    % Finally, put DIM back into DA and DB.
    if(allow_rect == 1)
        da = [dim(1,1),dim(2,1)];
        db = [dim(1,2),dim(2,2)];
    else
        da = dim(1,1);
        db = dim(1,2);
    end
end

% This function lets the user enter just a single number or a 1-by-2 vector
% for the dimensions, and then expands that number or vector to a full
% 2-by-2 matrix of dimensions appropriately.
function dim_mat = expand_dim(dim)
    sz = size(dim);
    if(max(sz) == 1) % user just entered a single number for DIM
        dim_mat = [dim,dim;dim,dim];
    elseif(min(sz) == 1) % user entered a 2-dimensional vector for DIM
        dim_mat = [dim(:).';dim(:).'];
    elseif(max(sz) == 2) % user entered a full 2-by-2 matrix for DIM
        dim_mat = dim;
    else
        error('expand_dim:InvalidDims','The dimensions must be provided in a matrix no larger than 2-by-2.');
    end
end