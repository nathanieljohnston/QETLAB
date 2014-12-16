%%  PAULICHANNEL   Generates a Pauli channel
%   This function has one required argument:
%     P: either a probability vector or a scalar indicating a number of
%        qubits
%
%   PHI = PauliChannel(P) is the Choi matrix of a Pauli channel. If P is a
%   scalar then P is a random P-qubit Pauli channel. If P is a probability
%   vector then it must have length 4^Q for some integer Q (which is the
%   number of qubits that PHI acts on), and the weight of the J-th Pauli
%   operator (in lexicographical order) in PHI's Kraus decomposition is
%   sqrt(P(J)). The output of this function is sparse.
%
%   URL: http://www.qetlab.com/PauliChannel

%   requires: ChoiMatrix.m, RandomProbabilities.m, update_odometer.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: December 16, 2014

function Phi = PauliChannel(p)

len_p = numel(p);
if(len_p == 1) % P is a scalar, so generate a random P-qubit Pauli channel
    q = p;
    len_p = 4^q;
    p = RandomProbabilities(len_p);
    
% Otherwise, do some error checking to make sure that P is a probability
% vector.
else
    % Let the user enter either a single vector P or a 4-by-4-by-...-by-4
    % array P.
    p = reshape(p,len_p,1);
    
    tol = len_p*eps^(3/4);
    if(isa(p,'cvx') == 0 && (min(p) < -tol || sum(p) > 1+tol))
        error('PauliChannel:InvalidP','P must be a probability vector: its entries must be non-negative and sum to 1.');
    end

    % Make sure that its length is a power of 2.
    [pow_of_2,double_q] = log2(len_p);
    if(abs(pow_of_2 - 1/2) > tol || mod(double_q,2) == 0)
        error('PauliChannel:InvalidP','P must contain 4^Q entries for some positive integer Q, since there are 4^Q Pauli operators on Q qubits.');
    end
    q = (double_q-1)/2; % number of qubits
end

% Finally, construct the channel.
Phi = sparse(len_p,len_p);
ind = zeros(1,q);
b4ones = 4*ones(1,q);
for j = 1:len_p
   Phi = Phi + p(j)*ChoiMatrix({Pauli(ind)});
   ind = update_odometer(ind,b4ones); % this loops over all base-4 strings of length Q
end