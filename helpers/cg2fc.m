function FC = CG2FC(CG,varargin)

%takes a Bell functional or behaviour CG that is
%in Collins-Gisin notation
% CG = [K   pB(0|1) p(0|2) ...
%	  pA(0|1) p(00|11) p(00|12) ...
%to Full Correlator notation
% FC =[  K   <B1> <B2> ... 
%      <A1> <A1B1> <A1B2> ...

%by default do the Bell functional, otherwise if behaviour==1 do the behaviour

if nargin == 1
	behaviour = 0;
else
	behaviour = varargin{1};
end

ia = size(CG,1)-1;
ib = size(CG,2)-1;

FC = zeros(ia+1,ib+1);

A = CG(2:end,1);
B = CG(1,2:end);
C = CG(2:end,2:end);

if behaviour == 0
	FC(1,1) = CG(1,1) + sum(A)/2 + sum(B)/2 + sum(sum(C))/4;
	FC(2:end,1) = A/2 + sum(C,2)/4;
	FC(1,2:end) = B/2 + sum(C,1)/4;
	FC(2:end,2:end) = C/4;
elseif behaviour == 1
	FC(1,1) = 1;
	FC(2:end,1) = 2*A-1;
	FC(1,2:end) = 2*B-1;
	FC(2:end,2:end) = (ones(ia,ib) - 2*A*ones(1,ib) - 2*ones(ia,1)*B + 4*C);
end
