function G = chshd(d)
%produces coefficients of the CHSH-d game

G = zeros(d,d,d,d);

for a = 0:d-1
	for b = 0:d-1
    	for x = 0:d-1
        	for y = 0:d-1
           		if(mod(a+b+x*y,d)==0)
             		G(a+1,b+1,x+1,y+1) = 1;
		        end
         	end
       	end
    end
end

G = G/d^2;
