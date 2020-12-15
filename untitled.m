clr;
%%
% H,W,T,D,R
% H = 10;
% W = 2;
% T = 4;
% D = 1;
% R = 1;
% 
% dD = 0.5*(W-D);
% dT = 0.5*(H-T);
% 
% R1 = sRectangle(0,W,0,H);
% R2 = sRectangle(0,dD,dT,dT+T);
% R3 = sRectangle(W-dD,W,dT,dT+T);
% C1 = sCircle(dD-R,dT,R);
% C2 = sCircle(dD-R,dT+T,R);
% C3 = sCircle(W-dD+R,dT,R);
% C4 = sCircle(W-dD+R,dT+T,R);
% 
% sdf = (R1 -(R2+R3) - (C1+C2+C3+C4)).';
% 
% sdf.show();
% axis off;

f1 = @(x) sqrt(x) - 0.5;
f2 = @(x) sqrt((x(:,1)+0.5).^2 + x(:,2).^2) - 0.5;

sdf1 = SDF(f1);
sdf2 = SDF(f2);

sdf1.BdBox = [-1,1,-1,1];
sdf2.BdBox = [-1.5,0.5,-1,1];

S = sdf1+sdf2;
S.show()

