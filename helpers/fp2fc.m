function FC = FP2FC(V,varargin)

%takes a Bell functional or behaviour V that is in 
%full probability notation V(a,b,x,y) to Full Correlator notation
% FC =[  K   <B1> <B2> ... 
%      <A1> <A1B1> <A1B2> ...
%where K is the independent term

%by default do the Bell functional, otherwise if behaviour==1 do the behaviour

if nargin == 1
	behaviour = 0;
else
	behaviour = varargin{1};
end

oa = size(V,1);
ob = size(V,2);
ia = size(V,3);
ib = size(V,4);

if oa ~= 2 || ob ~= 2
	disp('error');
	return
end

FC = zeros(1+ia,1+ib);

FC(1,1) = sum(V(:));

for x=1:ia
	FC(x+1,1) = sum(V(1,1,x,:) + V(1,2,x,:) - V(2,1,x,:) - V(2,2,x,:));
end

for y=1:ib
	FC(1,1+y) = sum(V(1,1,:,y) + V(2,1,:,y) - V(1,2,:,y) - V(2,2,:,y));
end

for x=1:ia
	for y=1:ib
		FC(x+1,y+1) = V(1,1,x,y) + V(2,2,x,y) - V(1,2,x,y) - V(2,1,x,y);
	end
end

if behaviour == 0

	FC = FC/4;
	
elseif behaviour == 1
	
	FC(1,1) = 1;
	FC(2:end,1) = FC(2:end,1)/ib;
	FC(1,2:end) = FC(1,2:end)/ia;
	
end

end
