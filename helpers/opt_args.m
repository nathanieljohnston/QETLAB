%%  OPT_ARGS    Handles optional input arguments for functions
%   This function has one required argument: DEF_ARGS, and an arbitrary
%   number of additional optional arguments
%
%   VARARGOUT = opt_args(DEF_ARGS,EXTRA_ARG1,EXTRA_ARG2,...) provides the
%   values EXTRA_ARG1, EXTRA_ARG2, ..., with the understanding that if a
%   particular EXTRA_ARGn is not provided as input, then DEF_ARGS(n) is
%   returned in its place. This is a helper function that is only used to
%   clean up and simplify code within other functions.
%
%   URL: http://www.qetlab.com/opt_args

%   requires: nothing
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: November 2, 2012

function varargout = opt_args(def_args,varargin)

% verify optional arguments
max_opt_args = length(def_args);
num_opt_args = length(varargin);
if num_opt_args > max_opt_args
    ST = dbstack(1);
    error(strcat(ST.name,':TooManyArguments'),'Received %d optional arguments, but expected no more than %d.',num_opt_args,max_opt_args);
end

% return argument values
varargout = def_args;
varargout(1:num_opt_args) = varargin;