clr;
%%
W = 90;
H = 350;

sdf = sRectangle(-W/2,W/2,-H,0);

msh = Mesh(sdf,'Quads',[4,30]);
msh = msh.generate();

msh.show();

%% 
fem = Fem(msh,'TimeStep',1/1e3,'TimeEnd',3);

fem.addSupport(fem.FindNodes('Top'),[1,1]);
fem.addGravity([0,-9810]);

%%
fem.Material    = NeoHookeanMaterial(1,0.33);
fem.Contraction = @(fem) Controller(fem);

fem.simulate();

%% movie
close all;
fem.set('Linestyle','none');
figure(101);

t = fem.Log.t;

for ii = 1:fps(t,60):numel(t)
    fem.set('Node',fem.Log.Node{ii});
    fem.show('Field',fem.Log.Stress{ii});
    axis(2*[-W W -H 5]);
    
    background();
    drawnow();
end

function dV = Controller(fem)
A   = -0.1;
Ctr = fem.get('Center0');
dV  = A*invlerp(Ctr(:,1),-15,15)*sin(30*fem.Time);
end