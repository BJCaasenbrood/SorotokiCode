clr;
%
f = sSphere(1);

S = f + f.translate([3,0,0]) + f.translate([3,0,0]);

f = S.smoothunion(S.',1);
f = f.smoothunion(S.'.',1);

f.BdBox =[-5,5,-5,5,-5,5];

obj = Gmodel(f,'Quality',20);
obj.bake.render();
view(30,30)