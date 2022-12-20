clr;
%% loading .stl file
obj = Gmodel('Bunny.stl');

%% rendering    
fig(101,[9,9]);
obj = obj.bake().render(); 
axis tight;

% material
obj.Texture = (egg);
view(0,10);
obj.update();

%
export_fig matcap_egg.png -q400 -a4