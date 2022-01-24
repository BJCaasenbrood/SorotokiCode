clr;
%% model
obj = Gmodel('Bunny.stl');
obj = obj.set('Shading','Face');
%% set texture
obj.bake().render(); 
obj = obj.fix;
%% showing section analysis
x = 80;
%
for ii = 1:80
    obj = obj.slice('Z',x);
    obj.update();
    view(30,30);
    x = x-1;
end