syms p11 p12 p13 p21 p22 p23 x

P1 = [p11;p12;p13];
P2 = [p21;p22;p23];

F = f(x,P1,P2);



diff(F,p11)







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F = f(x,P1,P2)
syms r1 r2
A = 5e-5;
r = 1;

H = [ 0,              0,             0;
      0, -0.5*r*sqrt(3), 0.5*r*sqrt(3);
     -r,          0.5*r,         0.5*r;
     -0,             -0,            -0;
      0,              0,             0;
      0,              0,             0];

h = 0.1;
    
F1 = H*P1(:)*rdelta(x-0.5,h);
F2 = H*P2(:)*(-rdelta(x-0.5,h)+rdelta(x-1+h,h));

%F1 = a1*H*P1(:)*(rdelta(x-0.5,h));
%F2 = a2*H*P2(:)*(-rdelta(x-0.5,h)+rdelta(x-1+h,h));

F = A*(F1 + F2);
end

function y = rdelta(x,h)
y = piecewise(abs(x) < h/2, 1/h, 1);
end
