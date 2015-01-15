%%  UPB  Generates an unextendible product basis
%   This function may be called in several different ways:
%
%   U = UPB(NAME) is a matrix containing as its columns the vectors in the
%   unextendible product basis specified by the string NAME. NAME must be
%   one of: 'GenShifts', 'Min4x4', 'Pyramid', 'QuadRes', 'Shifts',
%   'SixParam', or 'Tiles'. See the online documentation for descriptions
%   of these different UPBs.
%
%   [U,V,W,...] = UPB(NAME) is the same as above, except in this case, the
%   unextendible product basis is obtained by tensoring the columns of U,
%   V, W, ... together. That is, U, V, W, ... are the local vectors in the
%   unextendible product basis.
%
%   U = UPB(DIM) and [U,V,W,...] = UPB(DIM) are as above, except DIM is a
%   vector containing the local dimensions of the requested UPB rather than
%   the name of the UPB.
%
%   U = UPB(DIM,VERBOSE) and [U,V,W,...] = UPB(DIM,VERBOSE) are as above,
%   where VERBOSE is a flag (either 1 or 0, default 1) that indicates that
%   a reference for the returned UPB will or will not be displayed.
%
%   URL: http://www.qetlab.com/UPB

%   requires: FourierMatrix.m, IsTotallyNonsingular.m, MinUPBSize.m,
%             normalize_cols.m, one_factorization.m, opt_args.m,
%             opt_disp.m, perm_inv.m, PermuteSystems.m, sporth.m, Swap.m
%
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: November 12, 2014
%   URL: http://www.qetlab.com/UPB

function [u,varargout] = UPB(name,varargin)

    show_name = false;
    given_dims = false;
    revp = -1; % by default, don't permute systems around after we're done constructing the UPB
    if(isnumeric(name)) % user provided dimensions, not a name, so find an appropriate UPB
        % set optional argument defaults: verbose=1
        [verbose] = opt_args({ 1 },varargin{:});

        given_dims = true;
        np = length(name);
        
        % pre-load various references
        refs = {'K. Feng. Unextendible product bases and 1-factorization of complete graphs. Discrete Applied Mathematics, 154:942-949, 2006.', ...
                'D.P. DiVincenzo, T. Mor, P.W. Shor, J.A. Smolin, and B.M. Terhal. Unextendible product bases, uncompletable product bases and bound entanglement. Commun. Math. Phys. 238, 379-410, 2003.', ...
                'C.H. Bennett, D.P. DiVincenzo, T. Mor, P.W. Shor, J.A. Smolin, and B.M. Terhal. Unextendible product bases and bound entanglement. Phys. Rev. Lett. 82, 5385-5388, 1999.', ...
                'T.B. Pedersen. Characteristics of unextendible product bases. Thesis, Aarhus Universitet, Datalogisk Institut, 2002.', ...
                'N. Johnston. The minimum size of qubit unextendible product bases. In Proceedings of the 8th Conference on the Theory of Quantum Computation, Communication and Cryptography (TQC 2013). E-print: arXiv:1302.1604 [quant-ph], 2013.', ...
                'The realm of common sense (if there is only a single party, the only UPBs are full bases of the space).', ...
                'N. Alon and L. Lovasz. Unextendible product bases. J. Combinatorial Theory, Ser. A, 95:169-179, 2001.\nSee also: http://www.njohnston.ca/2013/03/how-to-construct-minimal-upbs/', ...
                'J. Chen and N. Johnston. The minimum size of unextendible product bases in the bipartite case (and some multipartite cases). Comm. Math. Phys., 333(1):351-365, 2015.'};

        [name,revp] = sort(name);
        if(np == 1)
            upbp{1} = eye(name);
            name = '';
            ref_ind = 6;
        elseif(np == 2 && min(name) <= 2)
            % In this case, there are no UPBs smaller than a basis of the
            % whole space.
            upbp{1} = repmat(eye(name(1)),1,name(2));
            upbp{2} = repmat(eye(name(2)),1,name(1));
            upbp{revp(1)} = Swap(repmat(eye(name(revp(1))),1,name(revp(2))).',[1,2],[name(revp(2)),name(revp(1))],1).';

            ref_ind = 2;
            name = '';
        elseif(np == 2 && all(name == [3,3]))
            name = 'Tiles';
            ref_ind = 3;
            show_name = true;
        elseif(np == 2 && all(name == [4,4]))
            name = 'Feng4x4';
            ref_ind = 1;
        elseif(np == 2 && all(name == name(1)) && mod(name(1),2) == 1 && isprime(2*name(1)-1))
            varargin = {name(1)};
            name = 'QuadRes';
            ref_ind = 2;
            show_name = true;
        elseif(np == 3 && all(name == [2,2,2]))
            name = 'Shifts';
            ref_ind = 3;
            show_name = true;
        elseif(np == 3 && all(name == [2,2,3]))
            name = 'Feng2x2x3';
            ref_ind = 1;
        elseif(np == 3 && all(name == [2,2,5]))
            name = 'Feng2x2x5';
            ref_ind = 1;
        elseif(np == 4 && all(name == [2,2,2,2]))
            name = 'Feng2x2x2x2';
            ref_ind = 1;
        elseif(np == 4 && all(name == [2,2,2,4]))
            name = 'Feng2x2x2x4';
            ref_ind = 1;
        elseif(np == 5 && all(name == [2,2,2,2,5]))
            name = 'Feng2x2x2x2x5';
            ref_ind = 1;
        elseif(np == 8 && all(name == [2,2,2,2,2,2,2,2]))
            name = 'John2^8';
            ref_ind = 5;
        elseif(mod(np,4) == 0 && all(name == 2*ones(1,np)))
            varargin = {np};
            name = 'John2^4k';
            ref_ind = 5;
        elseif(mod(np,2) == 1 && all(name == 2*ones(1,np)))
            varargin = {np};
            name = 'GenShifts';
            ref_ind = 2;
            show_name = true;
        elseif((mod(sum(name)-np,2) == 1 || sum(mod(name,2)==0) == 0) && sum(mod(name,2)==0) <= 1)
            varargin = {name};
            name = 'AlonLovasz';
            ref_ind = 7;
        elseif(name(end)-1 == sum(name(1:end-1)-1) && sum(name-1) >= 3 && mod(sum(name)-np,2) == 0)
            varargin = {name};
            name = 'CJBip';
            ref_ind = 8;
        elseif(np == 2 && all(name == [4,6]))
            name = 'CJBip46';
            ref_ind = 8;
        elseif(np == 3 && name(1) == 2 && name(2) == 2 && mod(name(3),4) == 1)
            varargin = {name(3)};
            name = 'CJ4k1';
            ref_ind = 8;
        else
            try
                min_size = MinUPBSize(name,0);
            catch err
                if(strcmpi(err.identifier,'MinUPBSize:MinSizeUnknown'))
                    error('UPB:MinSizeUnknown','No minimal UPB is currently known in the specified dimensions.');
                else
                    rethrow(err);
                end
            end
            error('UPB:HardToConstruct',['Minimal UPBs are known to have size ',num2str(min_size),' in this case, but their construction is complicated and not implemented by this script.']);
        end
    end

    if(strcmpi(name,'Shifts')) % GenShifts reduces to Shifts when there are 3 parties
        name = 'GenShifts';
        varargin = {3};
    end
    
    % the "Pyramid" UPB
    if(strcmpi(name,'Pyramid'))
        h = sqrt(1 + sqrt(5))/2;
        for j = 4:-1:0 % pre-allocate
            upbp{1}(:,j+1) = [cos(2*pi*j/5);sin(2*pi*j/5);h];
        end
        upbp{1} = 2*upbp{1}/sqrt(5+sqrt(5));
        upbp{2} = upbp{1}(:,[1,3,5,2,4]);

    % the "Tiles" UPB
    elseif(strcmpi(name,'Tiles'))
        upbp{1}(:,5) = ones(3,1)/sqrt(3); % pre-allocate
        upbp{1}(:,1) = [1;0;0];
        upbp{1}(:,2) = [1;-1;0]/sqrt(2);
        upbp{1}(:,3) = [0;0;1];
        upbp{1}(:,4) = [0;1;-1]/sqrt(2);
        upbp{2}(:,5) = ones(3,1)/sqrt(3);
        upbp{2}(:,1) = [1;-1;0]/sqrt(2);
        upbp{2}(:,2) = [0;0;1];
        upbp{2}(:,3) = [0;1;-1]/sqrt(2);
        upbp{2}(:,4) = [1;0;0];

    % the "Min4x4" UPB
    elseif(strcmpi(name,'Min4x4'))
        upbp{1}(:,8) = [0;0;1;0]; % pre-allocate
        upbp{1}(:,1) = [1;-3;1;1]/sqrt(12);
        upbp{1}(:,2) = [1;0;0;0];
        upbp{1}(:,3) = [0;1;2;1]/sqrt(6);
        upbp{1}(:,4) = [1;0;0;-1]/sqrt(2);
        upbp{1}(:,5) = [0;1;0;0];
        upbp{1}(:,6) = [3;1;-1;1]/sqrt(12);
        upbp{1}(:,7) = [0;1;1;0]/sqrt(2);
        upbp{2}(:,8) = [-1;1+sqrt(2);0;1]/sqrt(5+2*sqrt(2)); % pre-allocate
        upbp{2}(:,1) = [0;1;-3-sqrt(2);-1-sqrt(2)]/sqrt(15+8*sqrt(2));
        upbp{2}(:,2) = [1;0;0;0];
        upbp{2}(:,3) = [1;0;sqrt(2)-1;1]/sqrt(5-2*sqrt(2));
        upbp{2}(:,4) = [0;1;0;0];
        upbp{2}(:,5) = [-1;1+sqrt(2);0;1]/sqrt(5+2*sqrt(2));
        upbp{2}(:,6) = [0;0;1;0];
        upbp{2}(:,7) = [1;1;1;-sqrt(2)]/sqrt(5);

    % the "QuadRes" UPB
    elseif(strcmpi(name,'QuadRes'))
        if(isempty(varargin))
            error('UPB:InvalidArguments','When using NAME=QuadRes, you must specify a second input argument that gives the dimension of the desires UPB.')
        elseif(~isprime(2*varargin{1}-1))
            error('UPB:InvalidArguments','When using NAME=QuadRes, the second input argument DIM must be such that 2*DIM-1 is prime.')
        elseif(mod(varargin{1},2) == 0)
            error('UPB:InvalidArguments','When using NAME=QuadRes, the second input argument DIM must be odd.')
        end

        p = 2*varargin{1}-1;
        q = quad_residue(p);
        s = setdiff(1:p-1,q);
        sm = sum(exp(2i*pi*q/p));
        N = max(-sm,1+sm);
        
        % This UPB is obtained by choosing entries from the Fourier matrix
        % carefully.
        F = FourierMatrix(p);
        F(1,:) = sqrt(N)*F(1,:);
        upbp{1} = F([1,q+1],:);
        upbp{2} = F([1,mod(s(1)*q,p)+1],:);
        
        % normalize the output
        upbp{1} = upbp{1}./repmat(sqrt(sum(abs(upbp{1}).^2,1)),varargin{1},1);
        upbp{2} = upbp{2}./repmat(sqrt(sum(abs(upbp{2}).^2,1)),varargin{1},1);
    
    % the "SixParam" UPB (i.e., the one from Section IV.A of DMSST03)
    elseif(strcmpi(name,'SixParam'))
        if(isempty(varargin) || length(varargin{1}) ~= 6)
            error('UPB:InvalidArguments','When using NAME=SixParam, you must specify a vector containing the six parameters [gammaA,thetaA,phiA,gammaB,thetaB,phiB].')
        elseif(any(abs(sin(varargin{1}([1,2,4,5]))) < 10*eps) || any(abs(cos(varargin{1}([1,2,4,5]))) < 10*eps))
            error('UPB:InvalidArguments','When using NAME=SixParam, none of the gammaA, thetaA, gammaB, and thetaB parameters can be a multiple of pi/2.')
        end
        
        arg_cell = num2cell(varargin{1});
        [gammaA,thetaA,phiA,gammaB,thetaB,phiB] = arg_cell{:}; % give more memorable names to parameters

        NA = sqrt(cos(gammaA)^2 + sin(gammaA)^2 * cos(thetaA)^2);
        NB = sqrt(cos(gammaB)^2 + sin(gammaB)^2 * cos(thetaB)^2);
        
        upbp{1}(:,5) = [0;sin(gammaA)*cos(thetaA)*exp(1i*phiA);cos(gammaA)]/NA; % pre-allocate
        upbp{1}(:,1) = [1;0;0];
        upbp{1}(:,2) = [0;1;0];
        upbp{1}(:,3) = [cos(thetaA);0;sin(thetaA)];
        upbp{1}(:,4) = [sin(gammaA)*sin(thetaA);cos(gammaA)*exp(1i*phiA);-sin(gammaA)*cos(thetaA)];
        upbp{2}(:,5) = [0;sin(gammaB)*cos(thetaB)*exp(1i*phiB);cos(gammaB)]/NB;
        upbp{2}(:,1) = [0;1;0];
        upbp{2}(:,2) = [sin(gammaB)*sin(thetaB);cos(gammaB)*exp(1i*phiB);-sin(gammaB)*cos(thetaB)];
        upbp{2}(:,3) = [1;0;0];
        upbp{2}(:,4) = [cos(thetaB);0;sin(thetaB)];

    % the "GenShifts" UPB
    elseif(strcmpi(name,'GenShifts'))
        if(isempty(varargin))
            error('UPB:InvalidArguments','When using NAME=GenShifts, you must specify a second input argument that gives the number of parties in the desired UPB.')
        elseif(mod(varargin{1},2) == 0)
            error('UPB:InvalidArguments','When using NAME=GenShifts, the second input argument P must be odd.')
        end
        
        k = (varargin{1}+1)/2;
        upbp{1} = [cos((0:(1/k):(2-1/k))*pi/2);sin((0:(1/k):(2-1/k))*pi/2)];
        upbp{1} = upbp{1}(:,[1,k+1:-1:2,k+2:end]);
        
        for j = varargin{1}-1:-1:1 % pre-allocate memory
            upbp{j+1} = upbp{1}(:,[1,circshift(2:(2*k),[0,j])]);
        end
        
    % the "Feng2x2x2x2" UPB
    elseif(strcmpi(name,'Feng2x2x2x2'))
        b1 = eye(2);
        b2 = [1 1;1 -1]/sqrt(2);
        b3 = [cos(pi/3) sin(pi/3);-sin(pi/3) cos(pi/3)];
        
        upbp{1} = [b1(:,1),b1(:,2),b1(:,1),b2(:,1),b2(:,2),b2(:,1)];
        upbp{2} = [b1(:,1),b2(:,1),b1(:,2),b1(:,2),b2(:,2),b1(:,1)];
        upbp{3} = [b1(:,1),b2(:,1),b3(:,1),b2(:,2),b3(:,2),b1(:,2)];
        upbp{4} = [b1(:,1),b2(:,1),b3(:,1),b3(:,2),b1(:,2),b2(:,2)];

    % the "John2^8" UPB
    elseif(strcmpi(name,'John2^8'))
        b1 = eye(2);
        b2 = [1 1;1 -1]/sqrt(2);
        b3 = [cos(pi/3) sin(pi/3);-sin(pi/3) cos(pi/3)];
        b4 = [cos(pi/5) sin(pi/5);-sin(pi/5) cos(pi/5)];
        b5 = [cos(pi/7) sin(pi/7);-sin(pi/7) cos(pi/7)];
        
        upbp{1} = [b1(:,1),b2(:,1),b2(:,2),b1(:,2),b3(:,1),b4(:,1),b4(:,1),b3(:,1),b1(:,2),b3(:,2),b4(:,2)];
        upbp{2} = [b1(:,1),b1(:,2),b2(:,1),b3(:,1),b4(:,1),b4(:,1),b4(:,2),b4(:,2),b3(:,1),b2(:,2),b3(:,2)];
        upbp{3} = [b1(:,1),b2(:,1),b3(:,1),b3(:,2),b1(:,2),b3(:,2),b1(:,2),b4(:,1),b3(:,1),b2(:,2),b4(:,2)];
        upbp{4} = [b1(:,1),b2(:,1),b3(:,1),b3(:,1),b3(:,2),b4(:,1),b3(:,2),b2(:,2),b4(:,1),b4(:,2),b1(:,2)];
        upbp{5} = [b1(:,1),b2(:,1),b1(:,2),b3(:,1),b4(:,1),b2(:,2),b5(:,1),b3(:,2),b3(:,1),b5(:,2),b4(:,2)];
        upbp{6} = [b1(:,1),b2(:,1),b3(:,1),b2(:,2),b4(:,1),b4(:,2),b5(:,1),b5(:,2),b2(:,2),b1(:,2),b3(:,2)];
        upbp{7} = [b1(:,1),b2(:,1),b3(:,1),b4(:,1),b2(:,2),b4(:,2),b2(:,2),b1(:,2),b3(:,2),b5(:,1),b5(:,2)];
        upbp{8} = [b1(:,1),b2(:,1),b3(:,1),b4(:,1),b5(:,1),b1(:,2),b5(:,1),b3(:,2),b5(:,2),b4(:,2),b2(:,2)];

    % the "John2^4k" UPB
    elseif(strcmpi(name,'John2^4k'))
        if(isempty(varargin))
            error('UPB:InvalidArguments','When using NAME=John2^4k, you must specify a second input argument that gives the number of parties in the desired UPB.')
        elseif(mod(varargin{1},4) ~= 0 || varargin{1} <= 7)
            error('UPB:InvalidArguments','When using NAME=John2^4k, the second input argument P must equal 0 (mod 4) and must be at least 8.')
        end
        
        % Construct varargin{1}/2+2 distinct orthonormal bases of C^2.
        for j = (varargin{1}/2+2):-1:1 % pre-allocate
            jb = 2*pi/(2*j-1);
            b(:,:,j) = [cos(jb) sin(jb);-sin(jb) cos(jb)];
        end

        % the first 3 parties are special, so we do them separately
        upbp{1} = [reshape(repmat(reshape(b(:,1,1:varargin{1}/4+1),2,varargin{1}/4+1),2,1),2,varargin{1}/2+2), reshape(repmat(reshape(b(:,2,1:varargin{1}/4+1),2,varargin{1}/4+1),2,1),2,varargin{1}/2+2)];
        upbp{2} = [upbp{1}(:,1:varargin{1}/2+2),circshift(upbp{1}(:,varargin{1}/2+3:end),[0,2])];
        upbp{3} = [upbp{1}(:,1:varargin{1}/2+2),circshift(upbp{1}(:,varargin{1}/2+3:end),[0,4])];
        
        % now do the next 2k-4 parties
        if(varargin{1} > 8)
            upbp{4} = [reshape(b(:,1,:),2,varargin{1}/2+2), circshift(reshape(b(:,2,:),2,varargin{1}/2+2),[0,6])];
            upbp{5} = upbp{4}(:,[1:varargin{1}/2+2,reshape(fliplr(reshape(varargin{1}/2+3:varargin{1}+4,2,varargin{1}/4+1).').',1,varargin{1}/2+2)]);
            
            for j = 3:(varargin{1}/4-1)
                upbp{2*j} = [upbp{4}(:,1:varargin{1}/2+2),circshift(upbp{4}(:,varargin{1}/2+3:end),[0,2*(j-1)])];
                upbp{2*j+1} = [upbp{5}(:,1:varargin{1}/2+2),circshift(upbp{5}(:,varargin{1}/2+3:end),[0,2*(j-1)])];
            end
        end
        
        % Finally, do the last 2k+1 parties, which arise from finding a
        % 1-factorization of the complete graph.
        fac = one_factorization(varargin{1}/2+2);
        resb = reshape(b,2,varargin{1}+4);
        for j = 1:(varargin{1}/2+1)
            upbp{varargin{1}/2+j-1} = resb(:,[fac(j,:),varargin{1}/2+2+fac(j,:)]);
        end
        
    % the UPB from Theorem 3 of reference [8]
    elseif(strcmpi(name,'CJ4k1'))
        % Construct (varargin{1}+1)/2+1 distinct orthonormal bases of C^2.
        for j = ((varargin{1}+1)/2+1):-1:1 % pre-allocate
            jb = 2*pi/(2*j-1);
            b(:,:,j) = [cos(jb) sin(jb);-sin(jb) cos(jb)];
        end

        % the 2-dimensional parties are the same as in the Feng4m2 UPB
        upbp{1} = repmat(reshape(b(:,:,1:((varargin{1}+3)/4)),2,(varargin{1}+1)/2+1),1,2);
        u_ind = reshape((varargin{1}+1)/2+2:varargin{1}+3,2,(varargin{1}+3)/4);
        upbp{1}(:,(varargin{1}+1)/2+2:varargin{1}+3) = upbp{1}(:,reshape(u_ind([2,1],:),1,(varargin{1}+1)/2+1));
        b = reshape(b,2,varargin{1}+3);
        upbp{3} = CJ_Lemma6((varargin{1}-1)/4);
        upbp{2} = b(:,circshift(1:varargin{1}+3,[0,1]));
        upbp{2}(:,[(varargin{1}+3)/2,varargin{1}+3]) = upbp{2}(:,[varargin{1}+3,(varargin{1}+3)/2]);        

    % the UPB from Theorem 1 of reference [8]
    elseif(strcmpi(name,'CJBip'))
        maxD = varargin{1}(end);
        U = HollowUnitary(maxD);
        upbp{np} = [eye(maxD),U];

        if(maxD >= 20 && ~IsTotallyNonsingular(U,2:maxD-2))
            error('UPB:HardToConstruct',['Minimal UPBs are known to have size ',num2str(maxD*2),' in this case, but their construction is complicated and not implemented by this script.']);
        end
        
        prevdims = 0;
        for j = 1:np-1
            upbp{j} = CJ_Lemma5(maxD,varargin{1}(j)-1,1+prevdims);
            prevdims = prevdims + varargin{1}(j)-1;
        end

    % another UPB from Theorem 1 of reference [8]
    elseif(strcmpi(name,'CJBip46'))
        u = 1.64451358502312496885542269243;
        upbp{1} = normalize_cols([3 3 3 1 1 -1 -1 2 -2 0;2 1 1 2 3 -2 0 0 -2 -1;2 1 1 1 1 5 -3 -4 2 1;2 2 3 2 2 0 2 1 3 0]);
        upbp{2} = eye(6);
        upbp{2} = [upbp{2}(:,1:5),normalize_cols([0 0 (u-2)/(1-u) 1 (u-1)/(u-2);u*(u-2)/(2*u-3) 0 0 (3-2*u)/(u-2) 1;1 -1-u 0 0 1/(1+u);1 1 -u 0 0;0 u-1 1 1/(1-u) 0;u 1 1 1 1])];
        
    % the "Feng2x2x3" UPB
    elseif(strcmpi(name,'Feng2x2x3'))
        b1 = eye(2);
        b2 = [1 1;1 -1]/sqrt(2);

        upbp{3} = normalize_cols([1,1,2,1,0,0;0,0,-3,1,1,1;0,1,-1,-1,-3,0]);
        upbp{1} = [b1(:,1),b1(:,2),b1(:,1),b2(:,1),b2(:,2),b2(:,1)];
        upbp{2} = [b1(:,1),b2(:,1),b1(:,2),b1(:,2),b2(:,2),b1(:,1)];

    % the "Feng2x2x5" UPB
    elseif(strcmpi(name,'Feng2x2x5'))
        w = exp(pi*2i/3);
        b1 = eye(2);
        b2 = [1 1;1 -1]/sqrt(2);
        b3 = [cos(pi/3) sin(pi/3);-sin(pi/3) cos(pi/3)];
        b4 = [cos(pi/5) sin(pi/5);-sin(pi/5) cos(pi/5)];
        
        upbp{3} = normalize_cols([1,conj(w),w,0,0,0,1,0;0,w,conj(w),1,0,1,0,0;0,conj(w),0,0,0,1,w,1;0,0,w,0,1,1,conj(w),0;0,1,1,0,0,1,1,0]);
        upbp{1} = [b2(:,2),b2(:,1),b1(:,2),b1(:,1),b1(:,1),b1(:,2),b2(:,1),b2(:,2)];
        upbp{2} = [b4(:,2),b1(:,2),b4(:,1),b1(:,1),b3(:,1),b2(:,1),b3(:,2),b2(:,2)];

    % the "Feng4m2" UPB of Theorem 3.2 from the Feng paper
    % I don't know how to compute a 1-factorization of the complement graph
    % in this UPB's construction yet. The commented section below shows the
    % part of the construction that I *do* know how to implement.
%    elseif(strcmpi(name,'Feng4m2'))
%        if(isempty(varargin))
%            error('UPB:InvalidArguments','When using NAME=Feng4m2, you must specify a second input argument that gives the number of parties in the desired UPB.')
%        elseif(mod(varargin{1},4) ~= 2)
%            error('UPB:InvalidArguments','When using NAME=Feng4m2, the second input argument P must equal 2 (mod 4).')
%        end
%        
%        % Construct varargin{1}/2+1 distinct orthonormal bases of C^2.
%        for j = (varargin{1}/2+1):-1:1 % pre-allocate
%            jb = 2*pi/(2*j-1);
%            b(:,:,j) = [cos(jb) sin(jb);-sin(jb) cos(jb)];
%        end
%        
%        upbp{1} = repmat(reshape(b(:,:,1:((varargin{1}+2)/4)),2,varargin{1}/2+1),1,2); % these are the a_j's in the proof of the theorem
%        u_ind = reshape(varargin{1}/2+2:varargin{1}+2,2,(varargin{1}+2)/4);% we need to move the columns of upbp{1} around a little bit first
%        upbp{1}(:,varargin{1}/2+2:varargin{1}+2) = upbp{1}(:,reshape(u_ind([2,1],:),1,varargin{1}/2+1));
%        b = reshape(b,2,varargin{1}+2); % these are the vectors that will be used in all subsystems except the first
                
    % the "Feng4x4" UPB (Theorem 3.3(4) of Feng paper)
    elseif(strcmpi(name,'Feng4x4'))
        upbp{1}(:,8) = [1;-1;1;0]/sqrt(3); % pre-allocate
        upbp{1}(:,1:4) = eye(4);
        upbp{1}(:,5) = [0;1;1;1]/sqrt(3);
        upbp{1}(:,6) = [1;0;-1;1]/sqrt(3);
        upbp{1}(:,7) = [1;1;0;-1]/sqrt(3);
        upbp{2}(:,8) = [0;1;1;1]/sqrt(3);
        upbp{2}(:,[1,7,6,4]) = eye(4);
        upbp{2}(:,5) = [1;-1;1;0]/sqrt(3);
        upbp{2}(:,2) = [1;0;-1;1]/sqrt(3);
        upbp{2}(:,3) = [1;1;0;-1]/sqrt(3);
        
    % the "Feng2x2x2x4" UPB (Theorem 3.3(3) of Feng paper)
    elseif(strcmpi(name,'Feng2x2x2x4'))
        % need four different bases of C^2
        b1 = eye(2);
        b2 = [1 1;1 -1]/sqrt(2);
        b3 = [cos(pi/3) sin(pi/3);-sin(pi/3) cos(pi/3)];
        b4 = [cos(pi/5) sin(pi/5);-sin(pi/5) cos(pi/5)];
        
        upbp{4}(:,8) = [1;-1;1;0]/sqrt(3); % pre-allocate
        upbp{4}(:,1:4) = eye(4);
        upbp{4}(:,5) = [0;1;1;1]/sqrt(3);
        upbp{4}(:,6) = [1;0;-1;1]/sqrt(3);
        upbp{4}(:,7) = [1;1;0;-1]/sqrt(3);
        
        upbp{1} = [b1 b2 b3 b4];
        upbp{3} = [b1 b2 b3 b4];
        upbp{2} = [b1 b2 b3 b4];
        upbp{2} = upbp{2}(:,[1,3,7,5,4,2,6,8]);
        upbp{3} = upbp{3}(:,[1,5,7,3,4,8,6,2]);
        upbp{1} = upbp{1}(:,[1,3,7,5,8,6,2,4]);
        
    % the "Feng2x2x2x2x5" UPB (Theorem 3.3(5) of Feng paper)
    elseif(strcmpi(name,'Feng2x2x2x2x5'))
        w = exp(pi*2i/3);
        
        upbp{5}(:,10) = [1;conj(w);w;1;0]/2; % pre-allocate
        upbp{5}(:,1:5) = eye(5);
        upbp{5}(:,6) = [0;1;1;1;1]/2;
        upbp{5}(:,7) = [1;0;1;w;conj(w)]/2;
        upbp{5}(:,8) = [1;1;0;conj(w);w]/2;
        upbp{5}(:,9) = [1;w;conj(w);0;1]/2;

        % need five different bases of C^2
        upbp{1} = [eye(2), [1 1;1 -1]/sqrt(2), [cos(pi/3) sin(pi/3);-sin(pi/3) cos(pi/3)], [cos(pi/5) sin(pi/5);-sin(pi/5) cos(pi/5)], [cos(pi/7) sin(pi/7);-sin(pi/7) cos(pi/7)]];
        upbp{4} = upbp{1};
        upbp{3} = upbp{1};
        upbp{2} = upbp{1};
                
        upbp{2} = upbp{2}(:,[1,3,5,7,9,10,2,4,6,8]);
        upbp{3} = upbp{3}(:,[1,3,5,7,9,8,10,2,4,6]);
        upbp{4} = upbp{4}(:,[1,3,5,7,9,6,8,10,2,4]);
        upbp{1} = upbp{1}(:,[1,3,5,7,9,4,6,8,10,2]);
    elseif(strcmpi(name,'AlonLovasz'))
        dim = varargin{1};
        min_size = sum(dim) - np + 1;
        os = 0;
            
        if(mod(dim(1),2) == 0)
            upbp{1} = AL_even(dim(1),min_size);
        else
            upbp{1} = AL_odd(dim(1),min_size,os);
            os = os + (dim(1)-1)/2;
        end
        if(length(upbp{1}) == 1 && upbp{1} == -1)
            error('UPB:HardToConstruct',['Minimal UPBs are known to have size ',num2str(min_size),' in this case, but their construction is complicated and not implemented by this script.']);
        end

        for j = np:-1:2
            if(mod(dim(j),2) == 0)
                upbp{j} = AL_even(dim(j),min_size);                
            else
                upbp{j} = AL_odd(dim(j),min_size,os);
                os = os + (dim(j)-1)/2;
            end
            if(length(upbp{j}) == 1 && upbp{j} == -1)
                error('UPB:HardToConstruct',['Minimal UPBs are known to have size ',num2str(min_size),' in this case, but their construction is complicated and not implemented by this script.']);
            end
        end
    end
    
    if(length(revp) == 1)
        revp = 1:length(upbp);
    end
    revp = perm_inv(revp);
    u = upbp{revp(1)};
    varargout = upbp(revp(2:end));
    

    % If the user just requested one output argument, tensor local vectors
    % together.
    if(nargout <= 1 && length(varargout) >= 1)
        num_vec = size(u,2);
        % Do a clever column-by-column kronecker product that avoids having to
        % loop over every vector within each party. Note that this method
        % relies on not using sparse matrices.
        u = reshape(u,1,size(u,1),num_vec);
        for k = 1:length(varargout)-1
            v = reshape(varargout{k},size(varargout{k},1),1,num_vec);
            u = reshape(bsxfun(@times,v,u),1,size(u,2)*size(v,1),num_vec);
        end % we split the last iteration outside of the loop to reduce the number of required reshapes by 1
        v = reshape(varargout{end},size(varargout{end},1),1,num_vec);
        u = reshape(bsxfun(@times,u,v), size(u,2)*size(varargout{end},1),num_vec);
    end
    
    if(given_dims)
        if(show_name)
            opt_disp(['Generated the ',num2str(length(upbp{1})),'-state ''',name,''' UPB from:\n',refs{ref_ind},'\n'],verbose);
        else
            opt_disp(['Generated a minimal ',num2str(length(upbp{1})),'-state UPB from:\n',refs{ref_ind},'\n'],verbose);
        end
    end
end

% Computes all quadratic residues mod p (needed for the QuadRes UPB)
function q = quad_residue(p)
    q = unique(mod((1:floor(p/2)).^2,p));
end

% Computes the vectors in an Alon-Lovasz UPB in odd dimensions.
function W = AL_odd(r,c,os)

    its = 0;
    lin_indep = false;
    
    while ~lin_indep
        W = randn(r,c);
 
        for j = 2:c
            ind = intersect(1:j-1,mod([(j-(r-1)/2-os):(j-1-os),(j+1+os):(j+(r-1)/2+os)]-1,c)+1);
            tmp = null(W(:,ind)');
            W(:,j) = tmp*randn(size(tmp,2),1);
        end

        W = normalize_cols(W);
        lin_indep = IsTotallyNonsingular(W,r);
        its = its + 1;
        
        if(its >= 2 && ~lin_indep)
            W = -1;
            return
        end
    end
end

% Computes the vectors in an Alon-Lovasz UPB in one even dimension.
function W = AL_even(r,c)

    its = 0;
    lin_indep = false;
    
    while ~lin_indep
        r2 = c - r + 1;
        W = randn(r,c);

        for j = 2:c
            ind = setdiff(1:c,mod((j-(r2-1)/2:j+(r2-1)/2)-1,c)+1);
            tmp = null(W(:,ind)');
            W(:,j) = tmp*randn(size(tmp,2),1);
        end
        
        W = normalize_cols(W);
        lin_indep = IsTotallyNonsingular(W,r);
        its = its + 1;
        
        if(its >= 2 && ~lin_indep)
            W = -1;
            return
        end
    end
end

%%  CJ_LEMMA5    Construct a matrix W described by Lemma 5 of [CJ]
%   This function has three required argument:
%     Q: the positive integer q described by the lemma
%     R: the positive integer r described by the lemma
%     S: the positive integer s described by the lemma
%
%   W = CJ_Lemma5(Q,R,S) is a (R+1)-by-2Q matrix with columns of unit
%   length that satisfy the orthogonality condition (a) and the
%   nonsingularity condition (b) of Lemma 5.
%
%   This function has one optional argument:
%     NICE (default 0): a non-negative integer that specifies how "nice"
%                       the entries of W should be
%
%   W = CJ_Lemma5(Q,R,S,NICE) is the same as before if NICE = 0. If NICE is
%   positive integer, then W will be a symbolic matrix (rather than a
%   numeric matrix) whose entries are integers. In this case, the columns
%   of W are no longer scaled to have unit length, but rather are scaled to
%   make the entries the smallest integers possible. Lower values of NICE
%   lead to smaller entries, but in general take longer to compute (and may
%   not be able to be computed at all, if NICE is too small).
%
%   URL: http://www.njohnston.ca/publications/minimum-upbs/code/
%
%   References:
%   [CJ] J. Chen and N. Johnston. The Minimum Size of Unextendible Product
%        Bases in the Bipartite Case (and Some Multipartite Cases).
%        Preprint, 2012.

%   requires: IsTotallyNonsingular.m, normalize_cols.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   version: 1.00
%   last updated: December 19, 2012

function W = CJ_Lemma5(q,r,s,varargin)

    if(q < r + s)
        error('CJ_Lemma5:InvalidArgs','The arguments must satisfy q >= r + s.');
    end

    if(nargin == 3)
        nice = 0;
    else
        nice = varargin{1};
    end

    its = 0;
    satisfy_cond_b = false;
    while ~satisfy_cond_b
        % Construct a matrix W that satisfies condition (a) of Lemma 4.
        if(nice == 0)
            W = randn(r+1,2*q);
            W(:,q+1:(2*q)) = normalize_cols(W(:,q+1:(2*q)));
        else
            W = sym((1-2*floor(2*rand(r+1,2*q))).*(1+floor(nice*rand(r+1,2*q))));
        end
        lin_depen = 0;

        % Now generate the remaining conditions according to condition (b).
        for j = 0:q-1
            tmp_null = null(W(:,q+1+mod(j+(s:s+r-1),q)).');
            if(size(tmp_null,2) > 1)
                % Uh-oh, two of the randomly-generated columns are linearly
                % dependent: try again!
                lin_depen = 1;
                break;
            end
            W(:,j+1) = tmp_null;        

            % Get rid of the denominators in the new column, if the user wanted
            % all integers.
            if(nice > 0) 
                [~,cd] = numden(W(1,j+1));
                for k=2:size(W,1)
                    [~,den] = numden(W(k,j+1));
                    cd = lcm(cd,den);
                end

                W(:,j+1) = cd*W(:,j+1);
            end
        end

        % With probability 1, the matrix W also satisfies condition (b) of
        % Lemma 4. Nevertheless, we should make sure that this is the case, and
        % if it isn't, try again!
        if(lin_depen == 0)
            W = fliplr(W);
            satisfy_cond_b = IsTotallyNonsingular(double(W),r+1);
        end

        its = its + 1;
        if(its == 10 && nice ~= 0 && ~satisfy_cond_b)
            warning('CJ_Lemma5:LinearDependence', 'Tried (and failed) 10 times to find a matrix satisfying all imposed conditions. You may have better luck if you increase (or omit) the NICE argument.');
        end
    end
end

%%  CJ_LEMMA6    Construct a matrix W described by Lemma 6 of [CJ]
%   This function has one required argument:
%     K: the positive integer k described by the lemma
%
%   W = CJ_Lemma6(K) is a (4K+1)-by-(4K+4) matrix with columns of unit
%   length that satisfy the orthogonality conditions (i), (ii) and (iii)
%   and the nonsingularity condition (iv) of Lemma 6.
%
%   This function has one optional argument:
%     NICE (default 0): a non-negative integer that specifies how "nice"
%                       the entries of W should be
%
%   W = CJ_Lemma6(K,NICE) is the same as before if NICE = 0. If NICE is a
%   positive integer, then W will be a symbolic matrix (rather than a
%   numeric matrix) whose entries are integers. In this case, the columns
%   of W are no longer scaled to have unit length, but rather are scaled to
%   make the entries the smallest integers possible. Lower values of NICE
%   lead to smaller entries, but in general take longer to compute (and may
%   not be able to be computed at all, if NICE is too small).
%
%   URL: http://www.njohnston.ca/publications/minimum-upbs/code/
%
%   References:
%   [CJ] J. Chen and N. Johnston. The Minimum Size of Unextendible Product
%        Bases in the Bipartite Case (and Some Multipartite Cases).
%        Preprint, 2012.

%   requires: IsTotallyNonsingular.m, normalize_cols.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   version: 1.00
%   last updated: December 21, 2012

function W = CJ_Lemma6(k,varargin)

    d = 4*k + 1;

    if(nargin == 1)
        nice = 0;
    else
        nice = varargin{1};
    end

    its = 0;
    satisfy_cond_b = false;
    while ~satisfy_cond_b
        % Construct a matrix W that satisfies conditions (i), (ii) and (iii) of
        % Lemma 6.
        if(nice == 0)
            W = randn(d,d+3);
            W(:,1:2) = normalize_cols(W(:,1:2));
        else
            W = sym((1-2*floor(2*rand(d,d+3))).*(1+floor(nice*rand(d,d+3))));
        end

        for j = 2:d+2
            hlf = (d+3)/2;
            if(j < hlf)
                ind = setdiff(0:j-2,mod(j+1,hlf))+1;
            else
                ind = setdiff(0:j-1,[j-hlf,hlf+mod(j-1,hlf),hlf+mod(j+1,hlf)])+1;
            end
            tmp = null(W(:,ind).');
            W(:,j+1) = tmp(:,1);

            % Get rid of the denominators in the new column, if the user wanted
            % all integers.
            if(nice > 0) 
                [~,cd] = numden(W(1,j+1));
                for k=2:size(W,1)
                    [~,den] = numden(W(k,j+1));
                    cd = lcm(cd,den);
                end

                W(:,j+1) = cd*W(:,j+1);
            end
        end

        % With probability 1, the matrix W also satisfies condition (iv) of
        % Lemma 6. Nevertheless, we should make sure that this is the case, and
        % if it isn't, try again!
        satisfy_cond_b = IsTotallyNonsingular(double(W),d);

        its = its + 1;
        if(its == 10 && nice ~= 0 && ~satisfy_cond_b)
            warning('CJ_Lemma6:LinearDependence', 'Tried (and failed) 10 times to find a matrix satisfying all imposed conditions. You may have better luck if you increase (or omit) the NICE argument.');
        end
    end
end

%%  HOLLOWUNITARY    Produces a hollow unitary matrix
%   This function has one required argument:
%     DIM: a positive integer not equal to 1 or 3
%
%   U = HollowUnitary(DIM) is a unitary matrix with all of its diagonal
%   entries equal to 0 and all of its off-diagonal entries non-zero. Such a
%   matrix exists if and only if DIM = 2 or DIM >= 4.
%
%   URL: http://www.njohnston.ca/publications/minimum-upbs/code/
%
%   References:
%   [CJ] J. Chen and N. Johnston. The Minimum Size of Unextendible Product
%        Bases in the Bipartite Case (and Some Multipartite Cases).
%        Preprint, 2012.

%   requires: FourierMatrix.m
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   version: 1.00
%   last updated: December 19, 2012

function U = HollowUnitary(dim)

    if(dim == 3 || dim <= 1)
        error('HollowUnitary:InvalidDim','Hollow unitaries only exist when DIM = 2 or DIM >= 4.');
    end

    w2 = exp(2i*pi/(dim-2));
    F = FourierMatrix(dim-1);

    U2 = zeros(dim-1);
    for j = 1:dim-2
        U2 = U2 + w2^j*F(:,j+1)*F(:,j+1)';
    end

    U = [0,ones(1,dim-1)/sqrt(dim-1);ones(dim-1,1)/sqrt(dim-1),U2];
end