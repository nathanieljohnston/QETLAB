%%  XORGAMEVALUE    Computes the classical or quantum value of a nonlocal binary XOR game
%   This function has two required input arguments:
%     P: a matrix whose (s,t)-entry gives the probability that the referee
%        will give Alice the value s and Bob the value t
%     F: a binary matrix whose (s,t)-entry indicates the winning choice
%        (either 0 or 1) when Alice and Bob receive values s and t from the
%        referee
%
%   VAL = XORQuantumValue(P,F) is the classical value of the XOR game
%   specified by the probability matrix P and the binary matrix of winning
%   values F. That is, it is the optimal probability that Alice and Bob can
%   win the game if they are allowed determine a joint strategy beforehand,
%   but not allowed to communicate during the game itself.
%
%   This function has one optional input argument:
%     VTYPE (default 'classical'): one of 'classical' or 'quantum',
%     indicating what type of value of the game should be computed.
%
%   VAL = XORQuantumValue(P,F,VTYPE) is the value of the specified XOR game
%   in the setting where Alice and Bob can use strategies in the setting
%   (classical, quantum, or no-signalling) specified by VTYPE.
%
%   URL: http://www.qetlab.com/XORGameValue

%   requires: cvx (http://cvxr.com/cvx/)
%
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: March 6, 2015

function val = XORGameValue(p,f,varargin)

[s,t] = size(p);

% set optional argument defaults: VTYPE='classical'
[vtype] = opt_args({ 'classical' },varargin{:});

% do some error checking
tol = eps*s^2*t^2;
if(abs(sum(sum(p)) - 1) > tol)
    error('XORGameValue:InvalidP','P must be a probability matrix: its entries must sum to 1.');
elseif(min(min(p)) < -tol)
    error('XORGameValue:InvalidP','P must be a probability matrix: its entries must be nonnegative.');
elseif(~all([s,t] == size(f)))
    error('XORGameValue:InvalidDims','P and F must be matrices of the same size.');
end


% Compute the value of the game, depending on which type of value was
% requested.
if(strcmpi(vtype,'quantum'))
    % use semidefinite programming to compute the value of the game
    P = p.*(1-2*f);
    Q = [zeros(s),P;P',zeros(t)];

    cvx_begin sdp quiet
        cvx_precision default;
        variable X(s+t,s+t) symmetric
        maximize trace(Q*X)
        subject to
            diag(X) == 1;
            X >= 0;
    cvx_end

    % The above SDP actually computes the bias of the game. Convert it to the
    % value of the game now.
    val = real(cvx_optval)/4 + 1/2;
elseif(strcmpi(vtype,'classical'))
    val = 0; % at worst, our winning probability is 0... now try to improve

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

            val = max(val,p_win); % is this strategy better than other ones tried so far?
            if(val >= 1 - tol) % already optimal? quit
                return;
            end
        end
    end
else
    error('XORGameValue:InvalidVTYPE','VTYPE must be one of ''classical'' or ''quantum''.');
end