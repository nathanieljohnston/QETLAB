function V = FC2FP(FC,varargin)

%takes a Bell functional or behaviour FC that is in 
%Full Correlator notation
% FC =[  K   <B1> <B2> ... 
%      <A1> <A1B1> <A1B2> ...
%to full probability notation V(a,b,x,y)

%by default do the Bell functional, otherwise if behaviour==1 do the behaviour

%   author: Mateus Araújo
%   package: QETLAB
%   last updated: September 16, 2022

if nargin == 1
	behaviour = 0;
else
	behaviour = varargin{1};
end

ia = size(FC,1)-1;
ib = size(FC,2)-1;
V = zeros(2,2,ia,ib);

if behaviour == 0
	for x=1:ia
		for y=1:ib
			V(1,1,x,y) = FC(1,1)/(ia*ib) + FC(1+x,1)/ib + FC(1,1+y)/ia + FC(1+x,1+y);
			V(1,2,x,y) = FC(1,1)/(ia*ib) + FC(1+x,1)/ib - FC(1,1+y)/ia - FC(1+x,1+y);
			V(2,1,x,y) = FC(1,1)/(ia*ib) - FC(1+x,1)/ib + FC(1,1+y)/ia - FC(1+x,1+y);
			V(2,2,x,y) = FC(1,1)/(ia*ib) - FC(1+x,1)/ib - FC(1,1+y)/ia + FC(1+x,1+y);						
		end
	end
elseif behaviour == 1
	for x=1:ia
		for y=1:ib
			V(1,1,x,y) = 1 + FC(1+x,1) + FC(1,1+y) + FC(1+x,1+y);
			V(1,2,x,y) = 1 + FC(1+x,1) - FC(1,1+y) - FC(1+x,1+y);
			V(2,1,x,y) = 1 - FC(1+x,1) + FC(1,1+y) - FC(1+x,1+y);
			V(2,2,x,y) = 1 - FC(1+x,1) - FC(1,1+y) + FC(1+x,1+y);						
		end
	end
	V = V/4;
end

end
