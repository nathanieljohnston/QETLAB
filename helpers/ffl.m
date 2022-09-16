function G = ffl()
%produces coefficients of the Fortnow-Feige-Lovász game

p = [1 1; 1 0]/3;

G = zeros(2,2,2,2);
for a = 1:2
    for b = 1:2
        for x = 1:2
            for y = 1:2
                if (max(x,a) ~= max(y,b))
                    G(a,b,x,y) = p(x,y);
                 end
            end
        end
    end
end


end
