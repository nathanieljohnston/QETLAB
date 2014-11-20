%%  XORCLASSICALVALUE    Computes the classical value of a nonlocal binary XOR game
%   This function has two required input arguments:
%     P: a matrix whose (s,t)-entry gives the probability that the referee
%        will give Alice the value s and Bob the value t
%     F: a binary matrix whose (s,t)-entry indicates the winning choice
%        (either 0 or 1) when Alice and Bob receive values s and t from the
%        referee
%
%   WC = XORClassicalValue(P,F) is the classical value of the XOR game
%   specified by the probability matrix P and the binary matrix of winning
%   values F. That is, it is the optimal probability that Alice and Bob can
%   win the game if they are not allowed to communicate during the game.
%
%   URL: http://www.qetlab.com/XORClassicalValue

%   requires: nothing
%
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   version: 0.60
%   last updated: November 19, 2014

function wc = XORClassicalValue(p,f)

[s,t] = size(p);

% do some error checking
tol = eps*s^2*t^2;
if(abs(sum(sum(p)) - 1) > tol)
    error('XORClassicalValue:InvalidP','P must be a probability matrix: its entries must sum to 1.');
elseif(min(min(p)) < -tol)
    error('XORClassicalValue:InvalidP','P must be a probability matrix: its entries must be nonnegative.');
elseif(~all([s,t] == size(f)))
    error('XORClassicalValue:InvalidDims','P and F must be matrices of the same size.');
end
wc = 0; % at worst, our winning probability is 0... now try to improve

% Find the maximum probability of winning (this is NP-hard, so don't expect
% an easy way to do it: just loop over all strategies).
for a_ans = 0:2^s-1 % loop over Alice's answers
    for b_ans = 0:2^t-1 % loop over Bob's answers
        a_vec = bitget(a_ans,1:s);
        b_vec = bitget(b_ans,1:t);
        
        % Now compute the winning probability under this strategy: XOR
        % together Alice's responses and Bob's responses, then check where
        % the XORed value equals the value in the given matrix F. Where the
        % values match, multiply by the probability of getting that pair of
        % questions (i.e., multiply entry-wise by P) and then sum over the
        % rows and columns.
        p_win = sum(sum((mod(a_vec.'*ones(1,t) + ones(s,1)*b_vec,2) == f).*p));
        
        wc = max(wc,p_win); % is this strategy better than other ones tried so far?
        if(wc >= 1 - tol) % already optimal? quit
            return;
        end
    end
end