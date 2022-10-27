function CG = FP2CG(V,varargin)
%takes a Bell functional or behaviour V that is in 
%full probability notation V(a,b,x,y)
%to Collins-Gisin notation
% CG = [K   pB(0|1) pB(1|1) ... p(0|2) ...
%	  pA(0|1) p(00|11) p(01|11) ... p(00|12) ...

%by default do the Bell functional, otherwise if behaviour==1 do the behaviour

if nargin == 1
	behaviour = 0;
else
	behaviour = varargin{1};
end

[oa ob ia ib] = size(V);
alice_pars = ia*(oa-1)+1;
bob_pars = ib*(ob-1)+1;
aindex = @(a,x) (1 + a + (x-1)*(oa-1));
bindex = @(b,y) (1 + b + (y-1)*(ob-1));


CG = zeros(alice_pars,bob_pars);

if behaviour == 0

	CG(1,1) = sum(sum(V(oa,ob,:,:)));
	for a=1:oa-1
		for x=1:ia
			CG(aindex(a,x),1) = sum(V(a,ob,x,:) - V(oa,ob,x,:));
		end
	end
	for b=1:ob-1
		for y=1:ib
			CG(1,bindex(b,y)) = sum(V(oa,b,:,y) - V(oa,ob,:,y));
		end
	end
	for a=1:oa-1
		for b=1:ob-1
			for x=1:ia
				for y=1:ib
					CG(aindex(a,x),bindex(b,y)) = V(a,b,x,y) - V(a,ob,x,y) - V(oa,b,x,y) + V(oa,ob,x,y);
				end
			end
		end
	end

elseif behaviour == 1

	CG(1,1) = 1;

	for x=1:ia
		for a=1:oa-1
			CG(aindex(a,x),1) = sum(V(a,:,x,1));
		end
	end

	for y=1:ib
		for b=1:ob-1
			CG(1,bindex(b,y)) = sum(V(:,b,1,y));
		end
	end

	for x=1:ia
		for y=1:ib
			CG(aindex(1,x):aindex(oa-1,x),bindex(1,y):bindex(ob-1,y)) = V(1:oa-1,1:ob-1,x,y);
		end
	end

end


end




