clr;

R = 2;
W = 10;
H = 5;
T = 1.5;
D = .75;
N = 8;

S = sBellow(W,H,T,D,R);
S = S.repeat([0,H],N);
S = S.revolve();
S = S.translate([0,W,0]);
S = S.mirror([0,1,0]);
% 
S.BdBox = [-W,W,-25,25,-1,H*N+1];
% 
% obj = Gmodel(S,'Quality',50);
% obj.bake.render();
% 
% view(30,30);

S.show();