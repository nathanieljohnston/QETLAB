%%  UPDATE_ODOMETER    Increases a vector subject to limits on how large each entry can be
%
%   NEW_IND = update_odometer(OLD_IND,UPPER_LIM) increases the last entry
%   of the vector OLD_IND by 1, unless that would make it larger than the
%   last entry of the vector UPPER_LIM. In this case, it sets the last
%   entry to 0 and instead increases the *second*-last entry of OLD_IND,
%   unless that would make it larger than the second-last entry of
%   UPPER_LIM. In this case, it sets the second-last entry to 0 and instead
%   increases the *third*-last entry of OLD_IND (and so on; it works like
%   an odometer).
%
%   This function is useful when you want to have k nested for loops, but k
%   isn't specified beforehand. For example, instead for looping over i and
%   j going from 1 to 3, you could loop over a single variable going from
%   1 to 3^2, and set [i,j] = update_odometer([i,j],[3,3]) at each step
%   within the loop.
%
%   URL: http://www.qetlab.com/update_odometer

%   requires: nothing
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: December 11, 2014

function new_ind = update_odometer(old_ind,upper_lim)
    % Start by increasing the last index by 1.
    ind_len = length(old_ind);
    new_ind = old_ind;
    new_ind(end) = new_ind(end)+1;
        
    % Now we work the "odometer": repeatedly set each digit to 0 if it
    % is too high and carry the addition to the left until we hit a
    % digit that *isn't* too high.
    for j = ind_len:-1:1
        % If we've hit the upper limit in this entry, move onto the next
        % entry.
        if(new_ind(j) >= upper_lim(j))
            new_ind(j) = 0;
            if(j >= 2)
                new_ind(j-1) = new_ind(j-1) + 1;
            else % we're at the left end of the vector; just stop
                return;
            end
        else
            return; % always return if the odometer doesn't turn over
        end
    end
end