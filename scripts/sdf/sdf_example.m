clr;
%%
f = sSphere(1)/sCube(1.5/2);
% 
c = sCylinder(0,0,-1,1,0.5);
f = f - c - (c.') - (c.').';

obj = Gmodel(f,'Quality',80);
obj.bake.render();