%%  BELL    Produces a Bell state
%   This function has no required arguments.
%
%   PHI = Bell() is the Bell state (|0>|0> + |1>|1>)/sqrt(2)
%
%   This function has three optional arguments:
%     IND (default 0)
%     SP (default 0)
%     NRML (default 1)
%   
%   PHI = Bell(IND,SP,NRML) produces one of the following four Bell states,
%   depending on the value of IND:
%     0: (|0>|0> + |1>|1>)/sqrt(2)
%     1: (|0>|0> - |1>|1>)/sqrt(2)
%     2: (|0>|1> + |1>|0>)/sqrt(2)
%     3: (|0>|1> - |1>|0>)/sqrt(2)
%   The Bell state is sparse if SP = 1 and is full if SP = 0. It is
%   normalized to have Euclidean norm 1 if NRML = 1, and it is unnormalized
%   (i.e., each entry in the vector is 0 or 1) if NRML = 0.
%
%   URL: http://www.qetlab.com/Bell

%   requires: iden.m, opt_args.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: September 23, 2014

function phi = Bell(varargin)

% set optional argument defaults: ind=0, sp=0, nrml=1
[ind,sp,nrml] = opt_args({ 0, 0, 1 },varargin{:});

% construct an identity matrix, whose columns will be the local vectors
id = iden(2,sp);
colmod = mod(floor(ind/2),2); % we do (mod 2) so that values of ind > 3 work too

% construct the Bell state
phi = kron(id(:,1),id(:,1+colmod)) + ((-1)^ind)*kron(id(:,2),id(:,2-colmod));
if(nrml == 1)
    phi = phi/sqrt(2);
end