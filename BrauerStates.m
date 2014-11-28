%%  BRAUERSTATES    Produces all Brauer states
%   This function has two required arguments:
%     D: the local dimension
%     P: half of the number of parties (i.e., the states will live in
%     (2*P)-partite space)
%
%   B = BrauerStates(D,P) is a matrix whose columns are all of the
%   (unnormalized) "Brauer" states: states that are the P-fold tensor
%   product of the standard maximally-entangled state pure state on D local
%   dimensions. There are many such states, since there are many different
%   ways to group the 2*P parties into P pairs (with each pair
%   corresponding to one maximally-entangled state). The exact number of
%   such states is (2*P)!/(P!*2^P), which is the number of columns of B.
%
%   URL: http://www.qetlab.com/BrauerStates

%   requires: iden.m, MaxEntangled.m, opt_args.m, perfect_matchings.m,
%             PermuteSystems.m, Tensor.m
%
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: November 12, 2014

function B = BrauerStates(d,p)

phi = Tensor(MaxEntangled(d,1,0),p); % sparse, unnormalized

% The Brauer states are computed from perfect matchings of the complete
% graph. So compute all perfect matchings first.
pm = perfect_matchings(2*p);
npm = size(pm,1);
B = sparse(d^(2*p),npm);

% Now just turn these perfect matchings into the corresponding states.
for j = 1:npm
    B(:,j) = PermuteSystems(phi,pm(j,:),d*ones(1,2*p));
end