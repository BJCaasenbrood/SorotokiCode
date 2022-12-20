fig(101,[9,9]);
%% initial render
b = sSphere(1)/sCube(1.5/2);
b = b.render();
pause;
%% adding cylinders
c = sCylinder(0,0,-1,1,0.5);
a = c + c.rotate('x',90) +  c.rotate('y',90);
a = a.render('Texture',diffuse(0.1));
pause;
%% cutting cylinders from body
f = b - a;
cla; f = f.render('Quality',90);
