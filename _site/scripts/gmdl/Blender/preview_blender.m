clr;
%% preview
obj = Gmodel('Cube.stl');
obj.Texture = base;
%% set texture
figure(101); subplot(2,3,1);
obj.bake().render().update; 
obj.ground(); view(30,30); 
obj.update();

%% transform
Blender(obj,'Twist',{'z',30});
    subplot(2,3,3); obj.bake().render(); obj.ground();
    view(30,30); obj.update();

Blender(obj.reset(),'Scale',{'z',0.5}); 
    subplot(2,3,2); obj.bake().render(); obj.ground();
    view(30,30); obj.update();
    
Blender(obj.reset(),'Rotate',{'3D',30,20,30});
    subplot(2,3,4); obj.bake().render(); obj.ground();
    view(30,30); obj.update();

Blender(obj.reset(),'Translate',{'3D',0.5,0.5,0});
    subplot(2,3,5); obj.bake().render(); obj.ground();
    view(30,30); obj.update();
    
Blender(obj.reset(),'Curve',{'PCC',30,0});
    subplot(2,3,6); obj.bake().render(); obj.ground();    
    view(30,30); obj.update();
    