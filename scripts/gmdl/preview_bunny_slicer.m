clr;
%% model
obj = Gmodel('Bunny.stl');

%% set texture
obj.bake().render(); obj.ground;
background(gitpage);
view(20,20)

%% showing section analysis
gif('bunny_cut.gif','frame',gcf,'nodither');
x = 80;

for ii = 1:80
    obj.slice('Z',x);
    x = x-1;
    gif;
end