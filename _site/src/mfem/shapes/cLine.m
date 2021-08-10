
function c = cLine(x1,x2,y1,y2,z1,z2)   

c = @(P) [(1-P)*x1 + P*x2, (1-P)*y1 + P*y2, (1-P)*z1 + P*z2];


    