clc;
global R C L Rp u

R = .01;
C = 1e-3;
L = 1e-2;
Rp = 0.1;

u = -1;

[t,y]=ode15s(@ODE,[0 5],[0 0]);   

hold on
shade(t,y(:,2))

function dx = ODE(t,x)
global R C L Rp u
A = [-(1/R), -(1/C); 1/C, -Rp/(L*C)];
B = [1;0];
H = [x(1)/L; C*x(2)];

if t < 2,
    dx = A*H + B*(u);
else,
    dx = A*H;
end
   
end

