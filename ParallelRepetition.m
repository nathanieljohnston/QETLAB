function V2 = ParallelRepetition(V,rept)
%	produces the coefficients of rept parallel repetitions
%	of nonlocal game V given in full probability notation

%   requires: update_odometer.m
%
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: Sep 16, 2022

    oa = size(V,1);
    ob = size(V,2);
    ma = size(V,3);
    mb = size(V,4);
    
    V2 = V;
    
    % Are we playing the game multiple times? If so, adjust V
    % accordingly.
    if(rept > 1)

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
    end
    
end
