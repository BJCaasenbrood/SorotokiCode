function g = coadmap(x)
    W = [x(1);x(2);x(3)];
    U = x(4:6); 
    g = zeros(6);
    Wh = isomorph(W); 
    Uh = isomorph(U);
    
    g(1:3,1:3) = Wh;
    g(4:6,4:6) = Wh;
    g(1:3,4:6) = Uh;
end

function y = isomorph(x)
x1 = x(1); x2 = x(2); x3 = x(3);
y = [0, -x3, x2; x3, 0, -x1; -x2, x1, 0];
end
