%%  ISSEPARABLE    Determines whether or not a bipartite operator is separable
%   This function has one required argument:
%     X: a positive semidefinite matrix
%
%   SEP = IsSeparable(X) is either -1, 0, or 1. A value of 0 indicates that
%   X is entangled, a value of 1 indicates that X is separable, and a value
%   of -1 indicates that the separability of X could not be determined.
%
%   This function has four optional arguments:
%     DIM (default has both subsystems of equal dimension)
%     STR (default 2)
%     VERBOSE (default 1)
%     TOL (default eps^(3/8))
%
%   SEP = IsSeparable(X,DIM,STR,VERBOSE,TOL) is as above, where DIM is
%   a 1-by-2 vector containing the dimensions of the subsystems on which X
%   acts.
%
%   STR is an integer that determines how hard the script should work to
%   determine separability before giving up (STR = -1 means that the script
%   won't stop working until it finds an answer). Other valid values are 0,
%   1, 2, 3, ... In practice, if STR >= 4 then most computers will run out
%   of memory and/or the sun will explode before computation completes.
%
%   VERBOSE is a flag (either 0 or 1) that indicates that the script will
%   or will not display a line of text describing how it proved that X is
%   or is not separable.
%
%   TOL is the numerical tolerance used throughout the script.
%
%   URL: http://www.qetlab.com/IsSeparable

%   requires: ApplyMap.m, CVX (http://cvxr.com/cvx/), FilterNormalForm.m,
%             iden.m, InSeparableBall.m, IsPPT.m, IsPSD.m, jacobi_poly.m,
%             kpNorm.m, MaxEntangled.m, OperatorSchmidtDecomposition.m,
%             OperatorSchmidtRank.m, opt_args.m, opt_disp.m, PartialMap.m,
%             PartialTrace.m, PartialTranspose.m, PermutationOperator.m,
%             PermuteSystems.m, Realignment.m, SchmidtDecomposition.m,
%             SchmidtRank.m, sporth.m, Swap.m, SwapOperator.m,
%             SymmetricExtension.m, SymmetricInnerExtension.m,
%             SymmetricProjection.m, TraceNorm.m
%             
%   author: Nathaniel Johnston (nathaniel@njohnston.ca)
%   package: QETLAB
%   last updated: November 14, 2014

function sep = IsSeparable(X,varargin)

    X = full(X);
    if(~IsPSD(X))
        error('IsSeparable:NotPSD','X is not positive semidefinite, so the idea of it being separable does not make sense.');
    end
    
    lX = length(X);
    rX = rank(X);
    X = X/trace(X);
    sep = -1;

    % set optional argument defaults: dim=sqrt(length(X)), str=2, verbose=1, tol=eps^(1/4)
    [dim,str,verbose,tol] = opt_args({ round(sqrt(lX)), 2, 1, eps^(3/8) },varargin{:});
    if(str == -1)
        str = 1/eps; % keep going forever!
    end

    % allow the user to enter a single number for dim
    if(length(dim) == 1)
        dim = [dim,lX/dim];
        if abs(dim(2) - round(dim(2))) >= 2*lX*eps
            error('IsSeparable:InvalidDim','If DIM is a scalar, it must evenly divide length(X); please provide the DIM array containing the dimensions of the subsystems.');
        end
        dim(2) = round(dim(2));
    end
    nD = min(dim);
    xD = max(dim);
    pD = prod(dim);
    
    if(nD == 1)
        sep = 1;
        opt_disp(['Every positive semidefinite matrix is separable when one of the local dimensions is 1.\n'],verbose);
        return        
    end
    
    XA = PartialTrace(X,2,dim);
    XB = PartialTrace(X,1,dim);

    % pre-load various references
    refs = {'A. Peres. Separability criterion for density matrices. Phys. Rev. Lett., 77:1413-1415, 1996.', ... % refs[1]
            'M. Horodecki, P. Horodecki, and R. Horodecki. Separability of mixed states: Necessary and sufficient conditions. Phys. Lett. A, 223:1-8, 1996.', ...
            'P. Horodecki, M. Lewenstein, G. Vidal, and I. Cirac. Operational criterion and constructive checks for the separability of low-rank density matrices. Phys. Rev. A, 62:032310, 2000.', ...
            'K. Chen and L.-A. Wu. A matrix realignment method for recognizing entanglement. Quantum Inf. Comput., 3:193-202, 2003.', ...
            'F. Verstraete, J. Dehaene, and B. De Moor. Normal forms and entanglement measures for multipartite quantum states. Phys. Rev. A, 68:012103, 2003.', ...
            'K.-C. Ha and S.-H. Kye. Entanglement witnesses arising from exposed positive linear maps. Open Systems & Information Dynamics, 18:323-337, 2011.', ... % refs[6]
            'O. Gittsovich, O. Guehne, P. Hyllus, and J. Eisert. Unifying several separability conditions using the covariance matrix criterion. Phys. Rev. A, 78:052319, 2008.', ...
            'L. Gurvits and H. Barnum. Largest separable balls around the maximally mixed bipartite quantum state. Phys. Rev. A, 66:062311, 2002.', ...
            'H.-P. Breuer. Optimal entanglement criterion for mixed quantum states. Phys. Rev. Lett., 97:080501, 2006.', ...
            'W. Hall. Constructions of indecomposable positive maps based on a new criterion for indecomposability. E-print: arXiv:quant-ph/0607035, 2006.', ...
            'A. C. Doherty, P. A. Parrilo, and F. M. Spedalieri. A complete family of separability criteria. Phys. Rev. A, 69:022308, 2004.', ... % refs[11]
            'M. Navascues, M. Owari, and M. B. Plenio. A complete criterion for separability detection. Phys. Rev. Lett., 103:160404, 2009.', ...
            'N. Johnston. Separability from spectrum for qubit-qudit states. Phys. Rev. A, 88:062330, 2013.', ...
            'C.-J. Zhang, Y.-S. Zhang, S. Zhang, and G.-C. Guo. Entanglement detection beyond the cross-norm or realignment criterion. Phys. Rev. A, 77:060301(R), 2008.', ...
            'R. Hildebrand. Semidefinite descriptions of low-dimensional separable matrix cones. Linear Algebra Appl., 429:901-932, 2008,', ...
            'R. Hildebrand. Comparison of the PPT cone and the separable cone for 2-by-n systems. http://www-ljk.imag.fr/membres/Roland.Hildebrand/coreMPseminar2005_slides.pdf', ... % refs[16]
            'D. Cariello. Separability for weak irreducible matrices. E-print: arXiv:1311.7275 [quant-ph], 2013.', ...
            'L. Chen and D. Z. Djokovic. Separability problem for multipartite states of rank at most four. J. Phys. A: Math. Theor., 46:275304, 2013.', ...
            'G. Vidal and R. Tarrach. Robustness of entanglement. Phys. Rev. A, 59:141-155, 1999.'};

    % Start with the really easy separability checks (we always do these,
    % regardless of str).
    
    % PPT is used in so many tests, just always compute it
    ppt = IsPPT(X,2,dim,tol);
    
    % Check the PPT criterion
    if(~ppt)
        sep = 0;
        opt_disp(['Determined to be entangled via the PPT criterion. Reference:\n',refs{1},'\n'],verbose);
        return

    % Sometimes the PPT criterion is also sufficient for separability.
    elseif(pD <= 6 || min(dim) <= 1)
        sep = 1;
        opt_disp(['Determined to be separable via sufficiency of the PPT criterion in small dimensions. Reference:\n',refs{2},'\n'],verbose);
        return

    % Be careful with the next test! Checking that rX <= max(dim) is
    % *incorrect* here, and is not what is meant by refs{3}.
    elseif(rX <= 3 || rX <= rank(XB) || rX <= rank(XA))
        sep = 1;
        opt_disp(['Determined to be separable via sufficiency of the PPT criterion for low-rank operators. Reference:\n',refs{3},'\n'],verbose);
        return
    end
    
    % Realignment (aka computable cross norm) criterion. Only do it if
    % VERBOSE = 1, since we are about to do another test that is strictly
    % stronger.
    if(verbose && TraceNorm(Realignment(X,dim)) > 1)
        sep = 0;
        opt_disp(['Determined to be entangled via the realignment criterion. Reference:\n',refs{4},'\n'],verbose);
        return    
    end

    % Another test that is strictly stronger than the realignment criterion
    if(TraceNorm(Realignment(X - kron(XA,XB),dim)) > sqrt((1-trace(XA^2))*(1-trace(XB^2))))
        sep = 0;
        opt_disp(['Determined to be entangled by using Theorem 1 of reference:\n',refs{14},'\n'],verbose);
        return    
    end

    lam = sort(real(eig(X)),'descend'); % eigenvalues of X
    
    % There are some separability tests that work specifically in the
    % qubit-qudit (i.e., 2 \otimes n) case. Do these now.
    if(nD == 2)
        % check if X is separable from spectrum
        if((lam(1) - lam(2*xD-1))^2 <= 4*lam(2*xD-2)*lam(2*xD) + tol^2)
            sep = 1;
            opt_disp(['Determined to be separable by inspecting its eigenvalues. Reference:\n',refs{13},'\n'],verbose);
            return
        end
        
        % For the rest of the block matrix tests, we need the 2-dimensional
        % subsystem to be the *first* subsystem, so swap accordingly.
        if(dim(1) > 2)
            Xt = Swap(X,[1,2],dim);
        else
            Xt = X;
        end
        
        % Check if Lemma 1 of refs{13} applies to X. Also check the
        % Hildebrand 2xn results.
        A = Xt(1:xD,1:xD);
        B = Xt(1:xD,xD+1:2*xD);
        C = Xt(xD+1:2*xD,xD+1:2*xD);
        
        % This test is weaker than the next one, so only do it if VERBOSE =
        % 1.
        if(verbose && norm(B-B') < tol^2)
            sep = 1;
            opt_disp(['Determined to be separable by being a block Hankel matrix:\n',refs{15},'\n'],verbose);
            return
        end

        if(rank(B-B') <= 1 && ppt)
            sep = 1;
            opt_disp(['Determined to be separable by being a perturbed block Hankel matrix:\n',refs{16},'\n'],verbose);
            return
        end

        X_2n_ppt_check = [(5/6)*A-C/6, B; B', (5/6)*C-A/6];
        if(IsPSD(X_2n_ppt_check) && IsPPT(X_2n_ppt_check,2,[2,xD]))
            sep = 1;
            opt_disp(['Determined to be separable via the homothetic images approach of:\n',refs{16},'\n'],verbose);
            return
        end

        if(norm(B)^2 <= min(real(eig(A)))*min(real(eig(C))) + tol^2)
            sep = 1;
            opt_disp(['Determined to be separable by using Lemma 1 of reference:\n',refs{13},'\n'],verbose);
            return
        end
    end

    % There are conditions that are both necessary and sufficient when both
    % dimensions are 3 and the rank is 4
    if(rX == 4 && nD == 3 && xD == 3)
        % this method computes more determinants than are actually
        % necessary, but the speed loss isn't too great
        ranX = orth(X);
        for j = 6:-1:1
            for k = 7:-1:j+1
                for l = 8:-1:k+1
                    for m = 9:-1:l+1
                        p(j,k,l,m) = det(ranX([j,k,l,m],:));
                    end
                end
            end
        end
        
        F = det([p(1,2,4,5), p(1,3,4,6), p(2,3,5,6), p(1,2,4,6)+p(1,3,4,5), p(1,2,5,6)+p(2,3,4,5), p(1,3,5,6)+p(2,3,4,6);
                 p(1,2,7,8), p(1,3,7,9), p(2,3,8,9), p(1,2,7,9)+p(1,3,7,8), p(1,2,8,9)+p(2,3,7,8), p(1,3,8,9)+p(2,3,7,9);
                 p(4,5,7,8), p(4,6,7,9), p(5,6,8,9), p(4,5,7,9)+p(4,6,7,8), p(4,5,8,9)+p(5,6,7,8), p(4,6,8,9)+p(5,6,7,9);
                 p(1,2,4,8)-p(1,2,5,7), p(1,3,4,9)-p(1,3,6,7), p(2,3,5,9)-p(2,3,6,8), p(1,2,4,9)-p(1,2,6,7)+p(1,3,4,8)-p(1,3,5,7), p(1,2,5,9)-p(1,2,6,8)+p(2,3,4,8)-p(2,3,5,7), p(1,3,5,9)-p(1,3,6,8)+p(2,3,4,9)-p(2,3,6,7);
                 p(1,4,5,8)-p(2,4,5,7), p(1,4,6,9)-p(3,4,6,7), p(2,5,6,9)-p(3,5,6,8), p(1,4,5,9)-p(2,4,6,7)+p(1,4,6,8)-p(3,4,5,7), p(1,5,6,8)-p(2,5,6,7)+p(2,4,5,9)-p(3,4,5,8), p(1,5,6,9)-p(3,4,6,8)+p(2,4,6,9)-p(3,5,6,7);
                 p(1,5,7,8)-p(2,4,7,8), p(1,6,7,9)-p(3,4,7,9), p(2,6,8,9)-p(3,5,8,9), p(1,5,7,9)-p(2,4,7,9)+p(1,6,7,8)-p(3,4,7,8), p(1,5,8,9)-p(2,4,8,9)+p(2,6,7,8)-p(3,5,7,8), p(1,6,8,9)-p(3,4,8,9)+p(2,6,7,9)-p(3,5,7,9)]);
        
        sep = (abs(F) < max(tol^2,eps^(3/4))); % X is separable iff F is zero (or suffiently close to it, for numerical reasons)
        
        if(sep == 1)
            opt_disp(['Determined to be separable by checking the Chow form. Reference:\n',refs{18},'\n'],verbose);
        else
            opt_disp(['Determined to be entangled by checking the Chow form. Reference:\n',refs{18},'\n'],verbose);
        end
        return
    end
    
    % Check proximity of X with the maximally mixed state.
    if(InSeparableBall(X))
        sep = 1;
        opt_disp(['Determined to be separable by closeness to the maximally mixed state. Reference:\n',refs{8},'\n'],verbose);
        return    
    end
    
    % Check if X is a rank-1 perturbation of the identity, which is
    % necessarily separable if it's PPT, which we have already checked.
    if(lam(2) - lam(pD) < tol^2)
        sep = 1;
        opt_disp(['Determined to be separable by being a small rank-1 perturbation of the maximally-mixed state. Reference:\n',refs{19},'\n'],verbose);
        return   
    end
    
    % Check tensor rank equalling 2
    if(OperatorSchmidtRank(X,dim) <= 2)
        sep = 1;
        opt_disp(['Determined to be separable by having operator Schmidt rank at most 2. Reference:\n',refs{17},'\n'],verbose);
        return
    end

    % There is a family of known optimal positive maps in the qutrit-qutrit
    % case. Check for entanglement using these.
    if(dim(1) == 3 && dim(2) == 3)
        phi = MaxEntangled(3,0,0);
        for t=0:0.1:0.9
            for j=1:2
                if(t>0)
                    t = 1/t; % this is a weird way of using both t and 1/t as indices for the maps Phi we generate
                elseif(j>1)
                    break;
                end
                a = (1-t)^2/(1-t+t^2);
                b = t^2/(1-t+t^2);
                c = 1/(1-t+t^2);
                Phi = diag([a+1,c,b,b,a+1,c,c,b,a+1]) - phi*phi';

                if(~IsPSD(PartialMap(X,Phi,2,dim)))
                    sep = 0;
                    opt_disp(['Determined to be entangled via the positive map Phi[a,b,c] with a=',num2str(a),', b=',num2str(b),', c=',num2str(c),'. Reference:\n',refs{6},'\n'],verbose);
                    return
                end
            end
        end
    end

    % Use the Breuer-Hall positive maps (in even dimensions only) based on
    % antisymmetric unitary matrices.
    for p = 1:2
        if(mod(dim(p),2) == 0)
            phi = MaxEntangled(dim(p),0,0);
            U = kron(eye(dim(p)),fliplr(diag([ones(dim(p)/2,1);-ones(dim(p)/2,1)])));
            Phi = diag(ones(dim(p)^2,1)) - phi*phi' - U*SwapOperator(dim(p))*U';

            if(~IsPSD(PartialMap(X,Phi,p,dim)))
                sep = 0;
                opt_disp(['Determined to be entangled via the Breuer-Hall positive maps based on antisymmetric unitary matrices. References:\n',refs{9},'\n',refs{10},'\n'],verbose);
                return
            end
        end
    end
    
    % The next tests for entanglement and separability are slightly more
    % time-intensive, so we only do them if str >= 1.
    if(str >= 1)
        % Use the filter covariance matrix criterion (Filter CMC) for entanglement.
        % This strengthens the realignment criterion considerably, but we do the
        % realignment first anyway, since it is quicker and simpler (and thus
        % provides a more useful reference).
        try
            xi = FilterNormalForm(X,dim);

            if(sum(xi) > min( sqrt(nD*xD*(nD-1)*(xD-1)), nD*xD*(1 - 1/nD + (nD^2 - 1)/xD + min(0,(xD^2-nD^2)/xD - (xD-1)))/2 ) + 10*xD*eps)
                sep = 0;
                opt_disp(['Determined to be entangled via the Filter Covariance Matrix Criterion. Reference:\n',refs{7},'\n'],verbose);
                return    
            end

        % FilterNormalForm often is slow and results in an error when X is
        % low-rank, but the partial transpose criterion already took care of those
        % situations anyway. Nonetheless, we err on the side of caution and ignore
        % the cases when the state doesn't have a filter normal form.
        catch err
            if(~strcmpi(err.identifier,'FilterNormalForm:NoFNF'))
                rethrow(err);
            end
        end
        
        % Try the Guhne method of enhancing the separability within X.
        Xsep = X;
        XSeig = sort(real(eig(Xsep)));
        e_norm = norm(XSeig);
        it_ct = 0;
        while it_ct < 10
            lb = -1;
            % Iteratively find a product state whose maximal overlap with Xsep
            % is close to optimal. This has a non-zero chance of failing, so we
            % make sure that the obtained overlap isn't terrible before moving
            % on.
            while lb < XSeig(sum(dim)-1) - pD*eps;
                [lb,v] = sk_iterate(Xsep,1,dim);
            end

            % compute the minimum amount that we can subtract off v*v' and
            % still be positive semidefinite
            gen_eig = eig(Xsep,v*v');
            gen_eig = min(real(gen_eig(~isinf(gen_eig))));

            % actually subtract half that amount
            Xsep2 = Xsep - gen_eig*(v*v')/2;
            Xsep2 = Xsep2/trace(Xsep2);

            XSeig2 = sort(real(eig(Xsep)));
            e_norm2 = norm(XSeig);
            if(e_norm - e_norm2 < 10^(-4)) % don't use abs() -- we want e_norm2 to be smaller, not larger
                it_ct = it_ct + 1;
                if(e_norm - e_norm2 > -10^(-4)) % if e_norm2 isn't larger, then do the update anyway
                    XSeig = XSeig2;
                    Xsep = Xsep2;
                    e_norm = e_norm2;
                end
            else
                it_ct = 0;
                XSeig = XSeig2;
                Xsep = Xsep2;
                e_norm = e_norm2;

                if(InSeparableBall(Xsep))
                    sep = 1;
                    opt_disp('Determined to be separable by subtracting product states until the operator was close to the maximally-mixed state.',verbose);
                    return
                end

            end
        end
    end
    
    % If str is high enough, also try the stronger separability tests based on
    % symmetric extensions. These are much more time-intensive.
    for k = 2:str
        if(~SymmetricExtension(X,k,dim,1,1))
            sep = 0;
            opt_disp(['Determined to be entangled by not having a ',num2str(k),'-copy PPT symmetric extension. Reference:\n',refs{11},'\n'],verbose);
            return        
        end

        if(SymmetricInnerExtension(Xsep,k,dim,1)) % Use Xsep instead of X, since Xsep is in some sense "more" separable, so it is easier to detect
            sep = 1;
            opt_disp(['Determined to be separable via the semidefinite program based on ',num2str(k),'-copy PPT symmetric extensions from reference:\n',refs{12},'\n'],verbose);
            return
        end
    end
end