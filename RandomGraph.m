%%  RANDOMGRAPH    Generates a random graph
%   This function has one required argument:
%     N: the number of vertices in the graph
%
%   A = RandomGraph(N) generates the adjacency matrix of a random N-vertex
%       graph where each edge is present, independently of the others, with
%       probability 0.5.
%
%   This function has one optional argument:
%     P (default 0.5): a real number between 0 and 1; the probability of
%                      each edge being present
%
%   A = RandomGraph(N,P) generates the adjacency matrix of a random
%       N-vertex graph where each edge is present, independently of the
%       others, with probability P.
%
%   URL: http://www.qetlab.com/RandomGraph

%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: August 4, 2023

function A = RandomGraph(n,varargin)

    % set optional argument defaults: p=0.5
    [p] = opt_args({ 0.5 },varargin{:});

    if(p >= 1)% avoid division by 0 is upcoming formula
        A = ones(n) - eye(n);
    else
        A = triu(min(max(round(rand(n)/(2-2*p)),0),1),1);% random upper triangular part
        A = A + A';% symmetrize it
    end
end