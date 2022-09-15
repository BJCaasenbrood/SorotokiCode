clr;
%%
S = (sSphere(1).shell(0.1));
S = S.smoothdiff(sCube([-2,2,-2,2,0,2]),1);
obj = Gmodel(S,'Quality',120);
obj.bake.render();

