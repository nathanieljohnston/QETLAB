%%  RANDOMPOVM    Generates a random positive-operator valued measure
%   This function has two required arguments:
%     DIM: a scalar specifying the size of the POVM elements
%     OUTC: the number of outcomes of the measurement (i.e., the number of
%           POVM elements)
%
%   P = RandomPOVM(DIM,OUTC) generates a cell array with OUTC entries, each
%   of which is a DIM-by-DIM matrix that represents one of the measurement
%   operators in a POVM. Each element of the cell array is positive
%   semidefinite, and their sum is the identity matrix.
%
%   This function has one optional argument:
%     RE (default 0)
%
%   P = RandomPOVM(DIM,OUTC,RE) generates a cell array that represents a
%   POVM as above, and all of the entries of this POVM will be real if
%   RE=1.
%
%   URL: http://www.qetlab.com/RandomPOVM

%   requires: opt_args.m, RandomSuperoperator.m
%             
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: January 6, 2016

function P = RandomPOVM(dim,outc,varargin)

% set optional argument defaults: re=0
[re] = opt_args({ 0 },varargin{:});

% Generate the Choi matrix of a random quantum channel: we will generate
% the POVM from this.
Phi = RandomSuperoperator([outc,dim],0,1,re);

% The POVM elements are the diagonal blocks of the channel.
for j = outc:-1:1
    ind = (1:dim) + (j-1)*dim;
    P{j} = Phi(ind,ind);
end