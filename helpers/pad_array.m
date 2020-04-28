%%  Pad_array.m
%   Takes an array which will be prepadded with zeroes and an integer which
%   represents the number of zeroes to pad the array. This function takes three
%   arguments:
%     ARR: the array to be padded
%     PAD_NUM: the number of zeroes to prepad the array with
%     DIR: determines if the padding is before (0) or after (1) the input array
%
%   This function was written as a replacement for `padarray' within the
%   image processing package of MATLAB and Octave-Forge, as to not require having
%   those packages to use QETLAB.
%
%   URL: http://www.qetlab.com/pad_array

%   requires: nothing
%   author: Robert Materi
%   package: QETLAB
%   last updated: July 14, 2018

function bin = pad_array(arr, pad_num, dir)
    bin = arr;
    if (mod(pad_num,1) ~= 0 || pad_num < 0)
	    error ("pad_array: cannot have negative padding, expecting non-negative integer argument");
    end
    
    % prepend array
    if dir == 0 
        for d = 1:pad_num
            bin = [0, bin];
        end
    % postpend array
    elseif dir == 1
        for d = 1:pad_num
	        bin = [bin, 0];
        end
    else 
        error ("pad_array: end to pad array undefined, expecting 0 (prepend) or 1 (postpend)");
    end
end
