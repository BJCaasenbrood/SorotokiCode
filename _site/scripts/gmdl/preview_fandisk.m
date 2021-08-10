clr;
%% model
obj = Gmodel('Fandisk.stl');
obj.set('Shading','Face');
%% set texture
obj.Texture = diffuse(0.95);
obj.Emission = [0.15 0.75 0.15];

obj.bake();
    
%% set view
obj.render(); 
obj.update();
obj.ground();
axis tight;