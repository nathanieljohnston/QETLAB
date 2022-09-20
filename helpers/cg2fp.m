function V=CG2FP(CG,desc,varargin)

%takes a Bell functional or behaviour CG
%Collins-Gisin notation
% CG = [K   pB(0|1) p(0|2) ...
%	  pA(0|1) p(00|11) p(00|12) ...
%to full probability notation V(a,b,x,y)

%by default do the Bell functional, otherwise if behaviour==1 do the behaviour

%   author: Mateus Araújo
%   package: QETLAB
%   last updated: September 16, 2022

if nargin == 2
	behaviour = 0;
else
	behaviour = varargin{1};
end

oa = desc(1);
ob = desc(2);
ia = desc(3);
ib = desc(4);
aindex = @(a,x) (1 + a + (x-1)*(oa-1));
bindex = @(b,y) (1 + b + (y-1)*(ob-1));

if behaviour == 0

	for x=1:ia
		for y=1:ib
			for a=1:oa-1
				for b=1:ob-1
					V(a,b,x,y) = CG(1,1)/(ia*ib) + CG(aindex(a,x),1)/ib + CG(1,bindex(b,y))/ia + CG(aindex(a,x),bindex(b,y));
				end
			end
		end
	end
	
	for x=1:ia
		for y=1:ib
			for a=1:oa-1
				V(a,ob,x,y) = CG(1,1)/(ia*ib) + CG(aindex(a,x),1)/ib;
			end
		end
	end

	for x=1:ia
		for y=1:ib
			for b=1:ob-1
				V(oa,b,x,y) = CG(1,1)/(ia*ib) + CG(1,bindex(b,y))/ia;
			end
		end
	end
	
	for x=1:ia
		for y=1:ib
			V(oa,ob,x,y) =  CG(1,1)/(ia*ib);
		end
	end
	

elseif behaviour == 1

	for x=1:ia
		for y=1:ib
			V(1:oa-1,1:ob-1,x,y) = CG(aindex(1,x):aindex(oa-1,x),bindex(1,y):bindex(ob-1,y));
			V(1:oa-1,ob,x,y) = CG(aindex(1,x):aindex(oa-1,x),1) - sum(CG(aindex(1,x):aindex(oa-1,x),bindex(1,y):bindex(ob-1,y)),2);
			V(oa,1:ob-1,x,y) = CG(1,bindex(1,y):bindex(ob-1,y)) - sum(CG(aindex(1,x):aindex(oa-1,x),bindex(1,y):bindex(ob-1,y)),1);
			V(oa,ob,x,y) = CG(1,1) - sum(CG(aindex(1,x):aindex(oa-1,x),1),1) - sum(CG(1,bindex(1,y):bindex(ob-1,y)),2) + sum(sum(CG(aindex(1,x):aindex(oa-1,x),bindex(1,y):bindex(ob-1,y))));
		end
	end

end

end
