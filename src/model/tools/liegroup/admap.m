function g = admap(x)
    W = x(1:3);
    U = x(4:6); 
    g = zeros(6);
    Wh = isom(W); 
    Uh = isom(U);
    g(1:3,1:3) = Wh;
    g(4:6,4:6) = Wh;
    g(4:6,1:3) = Uh;
end
%--------------------------------------------------------------------------
function y = isom(v)
y = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0] ;
end