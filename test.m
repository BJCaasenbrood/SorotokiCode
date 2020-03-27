   v = [1,2,3,4,5,6];
tic
for ii = 1:1e6
   g = admap(v);
end
toc

%---------------------------------------- isomorhphism between SO(3) and R6
function y = isomSO3(x)
x1 = x(1); x2 = x(2); x3 = x(3);
y = [0, -x3, x2; x3, 0, -x1; -x2, x1, 0];
end

%----------------------------------------- adjoint action on lie alg. se(3)
function g = admap(x)
W = [x(1);x(2);x(3)];
U = x(4:6); 
g = zeros(6);
o = zeros(3);
Wh = isomSO3(W); 
Uh = isomSO3(U);
g = [Wh, o; Uh, Wh];
% g(1:3,1:3) = Wh;
% g(4:6,4:6) = Wh;
% g(4:6,1:3) = Uh;
end