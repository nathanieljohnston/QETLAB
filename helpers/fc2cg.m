function CG = FC2CG(FC,varargin)

%takes a Bell functional or behaviour FC that is in 
%Full Correlator notation
% FC =[  K   <B1> <B2> ... 
%      <A1> <A1B1> <A1B2> ...
%to Collins-Gisin notation
% CG = [K   pB(0|1) p(0|2) ...
%	  pA(0|1) p(00|11) p(00|12) ...

%by default do the Bell functional, otherwise if behaviour==1 do the behaviour

if nargin == 1
	behaviour = 0;
else
	behaviour = varargin{1};
end

ia = size(FC,1)-1;
ib = size(FC,2)-1;

CG = zeros(ia+1,ib+1);

A = FC(2:end,1);
B = FC(1,2:end);
C = FC(2:end,2:end);

if behaviour == 0
	CG(1,1) = FC(1,1) + sum(sum(C)) - sum(A) - sum(B);
	CG(2:end,1) = 2*A - 2*sum(C,2);
	CG(1,2:end) = 2*B - 2*sum(C,1);
	CG(2:end,2:end) = 4*C;
elseif behaviour == 1
	CG(1,1) = 1;
	CG(2:end,1) = (1+A)/2;
	CG(1,2:end) = (1+B)/2;
	CG(2:end,2:end) = (ones(ia,ib) + A*ones(1,ib) + ones(ia,1)*B + C)/4;
end
