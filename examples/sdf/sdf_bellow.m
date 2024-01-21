clr;

R1 = sRectangle([5,5],[20,15]);
R2 = sRectangle([0,0],[5,20]);
C1 = sCircle(5,[20,10]);

s = R1 + C1;
s = s.smoothunion(R2,10);

s = s.repeat([0,20],3);
S = s.revolve();
S = S.shell(1.5);
% S = S.translate([25,0,0]);
% S = S.mirror([1,0,0]);
% S.BdBox = [-50,50,-30,30,-5,65];
S = S.render('Quality',50);