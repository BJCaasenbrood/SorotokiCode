clr; 
%% Example: Sdf basic operations
% initial render
b = sSphere(1)/sCube(1.5/2);

% adding cylinders
c = sCylinder(0.5,[-1,1]);
a = c + c.rotate('x',90) +  c.rotate('y',90);
% a = a.render('Texture',matcap_diffuse(0.1));

% % cutting cylinders from body
f = b - a;
f = f.render('Quality',120);