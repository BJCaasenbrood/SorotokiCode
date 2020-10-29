clr;
%% model
obj = Gmodel('Bunny.stl');

%% set texture
obj.bake().render(); obj.ground; background(gitpage);
view(20,20)

%% showing section analysis
x = 80;

for ii = 1:80
    obj = obj.slice('Z',x);
    obj.update();
    x = x-1;
end