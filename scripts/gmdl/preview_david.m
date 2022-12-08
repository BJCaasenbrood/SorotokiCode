clr;
%% model
obj = Gmodel('David.stl');

%% set texture
obj.Texture = rot90(porcelain,2);

fig(101,[10,9]);
subplot(1,2,1);
patch('Vertices',obj.Node,'Faces',obj.Element,'EdgeColor','none');
lighting gouraud
material SHINY
view(210,10); axis equal;


subplot(1,2,2);
obj = obj.bake().render(); 
view(210,10); obj.update();

