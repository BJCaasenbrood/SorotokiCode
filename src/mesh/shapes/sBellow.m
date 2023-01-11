function sdf = sBellow(W,H,T,D,R)


eps = 1e-4*W;
sdf = Sdf(@(P) dBellow(P,W,H,T,D,R));
%sdf = Sdf(@(P) dBellow(P,R,W,H,T,D));
sdf.BdBox = [0,W,0,H];

% % generat sample points
% N   = 50;
% x   = linspace(-pi,pi,N).';
% S   = [(r-0.5*eps)*cos(x)+xc,...
%        (r-0.5*eps)*sin(x)+yc]; % sample set
%    
% sdf.Node    = S;
% sdf.Element = [(1:N-1).',(2:N).'];

end

function Dist = dBellow(P,W,H,T,D,r)
r2 = H/2 - D;
r1 = D + T;
x2 = min(max(W - r - r1 - r2,0),Inf);

C1 = dDiff(dCircle(P,r+r1,0,r1),dCircle(P,r+r1,0,r1-T));
C2 = dDiff(dCircle(P,r+r1+x2,H/2,r2),dCircle(P,r+r1+x2,H/2,r2-T));
C3 = dDiff(dCircle(P,r+r1,H,r1),dCircle(P,r+r1,H,r1-T));

R1 = dRectangle(P,r+r1,r+r1+x2,r1-T,r1);
R2 = dRectangle(P,r+r1,r+r1+x2,2*r2+D-T,2*r2+D);


Rc1 = dRectangle(P,r+r1,W,-D,D);
Rc2 = dRectangle(P,r+r1,W,-D+H,D+H);
Rc3 = dRectangle(P,0,r+r1+x2,r1,-D+H-T);

R = dUnion(R1,R2);
Rc = dUnion(dUnion(Rc1,Rc2),Rc3);
C = dUnion(dUnion(C1,C2),C3);

Dist = dIntersect(dDiff(dUnion(R,C),Rc),dRectangle(P,0,W,-1e-3,H+1e-3));
%Dist = dDiff(dUnion(R,C),Rc);
end


% 
% function Dist = dBellow(P,x0,x1,r1,r2,r3,r4,t)
% 
%   C1 = dCircle(P,r2+x0,0,r1);
%   C2 = dCircle(P,r2+x0,0,r2);
%   R1 = dRectangle(P,x0,r2+x0,0,r2);
%   R2 = dRectangle(P,r2+x0,r2+x1+x0,r1,r2);
%   C3 = dCircle(P,r2+x1+x0,r2+r3,r3);
%   C4 = dCircle(P,r2+x1+x0,r2+r3,r4);
%   R3 = dRectangle(P,r2+x1+x0,r2+x1+r4+x0,r1,r1+2*r4);
%   R4 = dRectangle(P,r2+x0,r2+x1+x0,r2+2*r3,r2+2*r3+t);
%   R6 = dRectangle(P,x0,r2+x0,r2+2*r3,r2+r2+2*r3);
%   C5 = dCircle(P,r2+x0,2*r2+2*r3,r2);
%   C6 = dCircle(P,r2+x0,2*r2+2*r3,r1);
%   
%   Dist0 = dDiff(dIntersect(R6,C5),C6);
%   Dist1 = dUnion(dDiff(dIntersect(R3,C4),C3),R4);
%   Dist  = dUnion(dUnion(dDiff(dUnion(R2,dIntersect(R1,...
%       C2)),C1),Dist1),Dist0);
% end