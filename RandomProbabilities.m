%%  RANDOMPROBABILITIES   Generates a random probability vector
%   This function has one required argument:
%     N: the length of the probability vector
%
%   P = RandomProbabilities(N) generates a length-N probability vector,
%   uniformly on the (N-1) unit simplex. That is, it generates a vector P
%   with N entries, each of which is between 0 and 1 and such that sum(P) =
%   1.
%
%   URL: http://www.qetlab.com/RandomProbabilities

%   requires: nothing
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: December 2, 2014

function p = RandomProbabilities(n)

% Start by splitting the interval [0,1] up into n bins and then find the
% width of each bin.
p = diff([0,sort(rand(1,n-1),'ascend'),1]);