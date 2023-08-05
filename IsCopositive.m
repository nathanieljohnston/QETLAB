%%  ISCOPOSITIVE    Determines whether or not a matrix is copositive
%   This function has one required input argument:
%     C: a matrix
%
%   ICP = IsCopositive(C) returns 0 (no), 1 (yes), or -1 (unsure),
%         indicating whether or not C is copositive. If an answer of -1 is
%         returned, use the optional arguments to tweak the function's
%         computational power.
%
%   This function has two optional input arguments:
%     K (optional, default 0): a non-negative integer that indicates the
%        level of the hierarchy used when testing copositivity
%     MODE (optional, default 'sos'): either 'sos' or 'nosdp', indicating
%        whether copositivity should be tested via the sum-of-squares
%        hierarchy or a no-SDP alternate hierarchy. The 'nosdp' hierarchy
%        is less accurate for a given level of the hierarchy, but is much
%        faster and less memory-intensive to run.
%
%   ICP = IsCopositive(C,K,MODE) returns 0 (no), 1 (yes), or -1 (unsure),
%         indicating whether or not C is copositive, by making use of the
%         K-th level of the hierarchy specified by MODE.
%
%   URL: http://www.qetlab.com/IsCopositive

%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: August 5, 2023

function icp = IsCopositive(C,varargin)

    % set optional argument defaults: K=0, MODE='sos'
    [k,mode] = opt_args({ 0, 'sos' },varargin{:});
    do_sos = strcmpi(mode,'sos');% true means SOS, false means SDP-free alternative

    n = length(C);% number of variables in the associated polynomial
    d = 2;% copositivity is determined by a degree-4 homogeneous polynomial, we use d = 4/2 = 2 is half the degree
    p = CopositivePolynomial(C);% the degree-4 polynomial associated with A

    if(do_sos)
        [lb,ub] = PolynomialSOS(p,n,d,k,'min',0);
    else
        [lb,ub] = PolynomialOptimize(p,n,d,k,'min',0);
    end

    if(lb >= -0.000000001)
        icp = 1;
    elseif(ub < 0)
        icp = 0;
    else
        icp = -1;
    end
end