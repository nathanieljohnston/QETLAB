%%  OPT_DISP  Display a message to the user (sometimes)
%   This function has two required arguments:
%     X: the object or message to be displayed
%     VERBOSE: a flag (1 or 0)
%
%   opt_disp(X,VERBOSE) calls fprintf(X) if VERBOSE = 1. If VERBOSE = 0
%   then this function does nothing.
%
%   Only use this function within other functions to easily display (or not
%   display) a message, if there is a "verbose" flag within that function.
%   If you always want to display a message, just use MATLAB's fprintf
%   function directly. The only purpose of this function is to improve
%   readability and reduce if statements in other code.
%
%   URL: http://www.qetlab.com/opt_disp

%   requires: nothing
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: December 3, 2012

function opt_disp(X,verbose)

if(verbose)
    fprintf(X);
end