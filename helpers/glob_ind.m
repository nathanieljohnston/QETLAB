%%  GLOB_IND    Creates a global index vector based on a local index vector
%   This function has two required arguments:
%     LI: a vector of local indices
%     DIM: the dimension of the local spaces (all of the same dimension, so
%          this is a scalar)
%
%   GI = glob_ind(LI,DIM) the global coordinate corresponding to the local
%   coordinates in LI. That is, GI is the position of the "1" in the vector
%   |LI(1)> \otimes |LI(2)> \otimes ...
%
%   URL: http://www.qetlab.com/glob_ind

%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: June 29, 2022

function gi = glob_ind(li,dim)
    gi = 1;
    num_li = length(li);
    for k = 1:num_li
        gi = gi + (dim^(k-1))*(li(k)-1);
    end
end