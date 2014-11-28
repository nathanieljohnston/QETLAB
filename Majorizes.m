%%  MAJORIZES    Determines whether or not a vector or matrix majorizes another
%   This function has two required arguments:
%     A: a vector or matrix
%     B: a vector or matrix
%
%   M = Majorizes(A,B) is either 1 or 0, indicating that A does or does not
%   majorize B. If A and B are matrices then majorization of their vectors
%   of singular values is checked.
%
%   This function has no optional arguments.
%
%   URL: http://www.qetlab.com/Majorizes

%   requires: nothing
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: March 4, 2014

function m = Majorizes(a,b)

if(min(size(a)) ~= 1) % if we received matrices as input, we check their singular values for majorization
    a = svd(full(a)); % automatically sorted in descending order
else
    a = sort(a(:),'descend');
end

% repeat for b
if(min(size(b)) ~= 1)
    b = svd(full(b));
else
    b = sort(b(:),'descend');
end

% pad with zeros if necessary
la = length(a);
lb = length(b);
if(la < lb)
    a = [a;zeros(lb-la,1)];
elseif(lb < la)
    b = [b;zeros(la-lb,1)];    
end

cta = 0;
ctb = -norm(a)*eps^(3/4); % for numerical robustness reasons

% now check majorization
for k = 1:length(a)
    cta = cta + a(k);
    ctb = ctb + b(k);
    if(cta < ctb)
        m = 0;
        return
    end
end

m = 1;