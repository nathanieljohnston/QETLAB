%%  ISUPB    Determines whether or not a set of product vectors form a UPB
%
%   IU = IsUPB(U,V,W,...) returns 1 or 0, indicating that the local states
%   specified by U, V, W, ... is or is not an unextendible product basis.
%   U, V, W, ... should be matrices, each with the same number of columns,
%   whose columns are the local vectors of the supposed UPB (i.e., each
%   matrix corresponds to one party).
%
%   This function has one optional output argument:
%     WIT
%
%   [IU,WIT] = IsUPB(U,V,W,...) is as above. If IU = 0 then WIT is a cell
%   containing local vectors for a product vector orthogonal to every
%   member of the non-UPB specified by U,V,W,.... If IU = 1 then WIT is
%   meaningless.
%
%   URL: http://www.qetlab.com/IsUPB

%   requires: opt_args.m, vec_partitions.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   last updated: October 21, 2014

function [iu,wit] = IsUPB(varargin)

wit = 0;
p = nargin; % number of parties in the UPB
s = size(varargin{1},2); % number of states in the UPB
for j = p:-1:1 % get local dimension of each party
    dim(j) = size(varargin{j},1);
end

pt = vec_partitions(1:s,p,dim-1); % all possible partitions of the states among the different parties
num_pt = length(pt); % number of partitions to loop through

if(num_pt == 0) % no partitions, which means that so few vectors were provided that it can't possibly be a UPB
    iu = 0;
    if(nargout > 1) % generate a witness, if requested
        orth_ct = 0;
        wit = cell(1,p);
        for l = p:-1:1
            if(orth_ct >= s) % the witness we have generated is already orthogonal to all local vectors: just pick the first basis state from here on
                wit{l} = zeros(dim(l),1);
                wit{l}(1) = 1;
            else
                ct_upd = min(dim(l)-1,s-orth_ct);
                wit{l} = null(varargin{l}(:,orth_ct+1:orth_ct+ct_upd).'); % these are the states orthogonal to the local states provided
                wit{l} = wit{l}(:,1); % we only need one witness state
                orth_ct = orth_ct + ct_upd;
            end
        end
    end
    return
end

% There are non-trivial partitions to loop through. Do it!
for k = 1:num_pt % loop through the partitions
    num_orth = 1;
    for l = p:-1:1 % loop through the parties
        V = varargin{l}(:,pt{k}{l}); % these are the local states that we are choosing on the l-th party
        orth_state{l} = null(V.'); % these are the states orthogonal to the local states
        num_orth = num_orth*size(orth_state{l},2);
        if(num_orth >= 1)
            orth_state{l} = orth_state{l}(:,1); % we just do this so that the wit variable is prettier
        end
    end
            
    if(num_orth >= 1) % all parties have codimension at least 1: not a UPB
        iu = 0;
        wit = orth_state;
        return
    end
end
    
iu = 1;