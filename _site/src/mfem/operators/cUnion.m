
function c = cUnion(c1,c2)   

c = @(x) 0.5*((1-sign(x-0.5)).*c1(2*x) + (1+sign(x-0.5)).*c2(2*x - 1));
end


    