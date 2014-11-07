%%  BRAUERSTATES    Produces all Brauer states
%   This function has two required arguments:
%     K: half of the number of parties (i.e., the states will live in
%     (2*K)-partite space)
%     N: the local dimension
%
%   B = BrauerStates(K,N) is a matrix whose columns are all of the
%   (unnormalized) "Brauer" states: states that are the K-fold tensor
%   product of the standard maximally-entangled state pure state on N local
%   dimensions. There are many such states, since there are many different
%   ways to group the 2*K parties into K pairs (with each pair
%   corresponding to one maximally-entangled state). The exact number of
%   such states is (2*K)!/(K!*2^K), which is the number of columns of B.
%
%   URL: http://www.qetlab.com/BrauerStates

%   requires: iden.m, MaxEntangled.m, opt_args.m, perfect_matchings.m,
%             PermuteSystems.m, Tensor.m
%
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   version: 0.50
%   last updated: November 6, 2014

function B = BrauerStates(k,n)

phi = Tensor(MaxEntangled(n,1,0),k); % sparse, unnormalized

% The Brauer states are computed from perfect matchings of the complete
% graph. So compute all perfect matchings first.
pm = perfect_matchings(2*k);
npm = size(pm,1);
B = sparse(n^(2*k),npm);

% Now just turn these perfect matchings into the corresponding states.
for j = 1:npm
    B(:,j) = PermuteSystems(phi,pm(j,:),n*ones(1,2*k));
end