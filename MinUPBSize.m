%%  MINUPBSIZE  Gives the minimum cardinality of an unextendible product basis in given dimensions
%   This function has one required input argument:
%     DIM: a vector containing the local dimensions
%
%   S = MinUPBSize(DIM) is the minimum possible number of vectors in an
%   unextendible product basis in a Hilbert space whose local dimensions
%   are given by DIM. A reference to a journal article that provides a
%   proof of the minimal size claimed by this script is also be displayed.
%
%   This function has one optional input argument:
%     VERBOSE (default 1)
%     
%   S = MinUPBSize(DIM,1) is as above. S = MinUPBSize(DIM,0) produces the
%   same output, but does not print the reference to a journal article.
%
%   URL: http://www.qetlab.com/MinUPBSize

%   requires: opt_args.m, opt_disp.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: June 22, 2013

function s = MinUPBSize(dim,varargin)

% set optional argument defaults: verbose=1
[verbose] = opt_args({ 1 },varargin{:});

% pre-load various references
refs = {'K. Feng. Unextendible product bases and 1-factorization of complete graphs. Discrete Applied Mathematics, 154:942–949, 2006.', ...
        'D.P. DiVincenzo, T. Mor, P.W. Shor, J.A. Smolin, and B.M. Terhal. Unextendible product bases, uncompletable product bases and bound entanglement. Commun. Math. Phys. 238:379–410, 2003.', ...
        'N. Alon and L. Lovasz. Unextendible product bases. J. Combinatorial Theory, Ser. A, 95:169-179, 2001.', ...
        'J. Chen and N. Johnston. The minimum size of unextendible product bases in the bipartite case (and some multipartite cases). Comm. Math. Phys., 333(1):351-365, 2015.', ...
        'N. Johnston. The minimum size of qubit unextendible product bases. In Proceedings of the 8th Conference on the Theory of Quantum Computation, Communication and Cryptography (TQC 2013). E-print: arXiv:1302.1604 [quant-ph], 2013.'};
        
if(verbose > 1 && all(size(dim) == [1,1]))
    error('MinUPBSize:IncorrectVerbose',['VERBOSE should be 0 or 1, not ',num2str(verbose),'. You perhaps meant to call this function like this: MinUPBSize([',num2str(dim),',',num2str(verbose),'])']);
end

np = length(dim);
dim = sort(dim);
tlb = sum(dim) - np + 1; % trivial lower bound on s
x_dim = max(dim);
    
if(np == 2 && min(dim) <= 2)
    ref_ind = 2;
    s = prod(dim);
elseif(mod(tlb,2) == 0 || all(mod(dim,2) == 1)) % if tlb is even or all dims are odd
    ref_ind = 3;
    s = tlb;
elseif(np == 2 && x_dim-1 >= tlb-x_dim && tlb-x_dim >= 3)
    ref_ind = 4;
    s = tlb + 1;
elseif((np == 4 && all(dim == [2,2,2,2])) || (mod(np,4) == 2 && all(dim == 2*ones(1,np))))
    ref_ind = 1;
    s = tlb + 1;
elseif(np == 8 && all(dim == [2,2,2,2,2,2,2,2]))
    ref_ind = 5;
    s = 11;
elseif(all(dim == 2*ones(1,np)))
    ref_ind = 5;
    s = tlb + 3;
elseif(np == 3 && (all(dim == [2,2,3]) || all(dim == [2,2,5])))
    ref_ind = 1;
    s = tlb + 1;
elseif(np == 3 && dim(1) == 2 && dim(2) == 2 && mod(dim(3),4) == 1)
    ref_ind = 4;
    s = tlb + 1;
else
    error('MinUPBSize:MinSizeUnknown','The minimum size of UPBs in the specified space is not currently known.');
end
    
opt_disp(['A proof of the minimal size in this case can be found in:\n',refs{ref_ind},'\n'],verbose);