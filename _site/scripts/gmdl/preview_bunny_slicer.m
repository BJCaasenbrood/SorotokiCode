clr;
%% model
obj = Gmodel('Soft_ModulePlanar_TPU.stl');
obj = obj.set('Shading','Face');
%% set texture
obj = Blender(obj,'Rotate',{'y',90});
obj.bake().render(); 
obj = obj.fix;

view(90,0);

%% showing section analysis
x = 80;
% 
% for ii = 1:80
     obj = obj.slice('Z',-6);
%     obj.update();
%     x = x-1;
% end