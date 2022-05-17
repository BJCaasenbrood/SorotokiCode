clr;
%% model
sdf = sCube(-1,1,-1,1,-1,1);
obj = Gmodel(sdf,'Quality',30,'Shading','Face');

%% set texture
obj.set('Texture',softmath);
obj.bake().render(); 
   