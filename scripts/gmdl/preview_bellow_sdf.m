clr;
%% model
Dist = @(x) SDF(x);
BdBox =[-50,80,-50,50,0,132];

obj = Gmodel(Dist,BdBox,'Quality',70);

%% set texture
obj.Texture = studioclay;
obj.bake();

%% show
obj.render(); 

function Dist = SDF(x)
x  = pRepeat(x,[0,0,22]);
y  = dRevolve(x);
B1 = Bellow(y,5,4,6,5,7,5,2);
y  = dTranslate(x,[44,0,0]);
y  = dRevolve(y);
B2 = Bellow(y,5,4,6,5,7,5,2);
Dist = dSmoothUnion(B1,B2,7);
end

function Dist = Bellow(P,r0,r1,r2,r3,r4,x,t)
  C1 = dCircle(P,r2+r0,0,r1);
  C2 = dCircle(P,r2+r0,0,r2);
  R1 = dRectangle(P,r0,r2+r0,0,r2);
  R2 = dRectangle(P,r2+r0,r2+x+r0,r1,r2);
  C3 = dCircle(P,r2+x+r0,r2+r3,r3);
  C4 = dCircle(P,r2+x+r0,r2+r3,r4);
  R3 = dRectangle(P,r2+x+r0,r2+x+r4+r0,r1,r1+2*r4);
  R4 = dRectangle(P,r2+r0,r2+x+r0,r2+2*r3,r2+2*r3+t);
  R6 = dRectangle(P,r0,r2+r0,r2+2*r3,r2+r2+2*r3);
  C5 = dCircle(P,r2+r0,2*r2+2*r3,r2);
  C6 = dCircle(P,r2+r0,2*r2+2*r3,r1);
  
  Dist0 = dDiff(dIntersect(R6,C5),C6);
  Dist1 = dUnion(dDiff(dIntersect(R3,C4),C3),R4);
  Dist  = dUnion(dUnion(dDiff(dUnion(R2,dIntersect(R1,C2)),C1),Dist1),Dist0);
end
