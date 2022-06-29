%%  Dec_To_Bin 	Takes a decimal number and turns it to a vector representing the 
%   binary version with the MSB to the right. This function has one required 
%   argument:
%     DEC: the number to convert to the binary form
%
%   This function was written as a replacement for `di2bi' within the
%   communications package of MatLab and Octave-Forge, as to not require having
%   those packages to use QETLAB.
%
%   URL: http://www.qetlab.com/dec_to_bin

%   requires: nothing
%   author: Robert Materi
%   package: QETLAB
%   last updated: June 25, 2018

function bin = dec_to_bin(dec)
    bin = [];
    if (mod(dec,1) ~= 0 || dec < 0)
	    error ("dec_to_bin: expecting non-negative integer argument");
    end
    while dec >= 1
	    bin = [bin, mod(dec,2)];
	    dec = floor(dec/2);
    end
end
