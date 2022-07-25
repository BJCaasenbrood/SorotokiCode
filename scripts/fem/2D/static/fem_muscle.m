clr;
%%
W = 30;
H = 350;

sdf = sRectangle(-W/2,W/2,-H,0);

msh = Mesh(sdf,'Quads',[2,60]);
msh = msh.generate();

msh.show();

%% 
fem = Fem(msh,'TimeStep',1/150,'TimeEnd',3);

fem.AddConstraint('Support',fem.FindNodes('Top'),[1,1]);
fem.AddConstraint('Gravity',[],[0,-9810]);

%%
fem.Material    = Dragonskin10(100);
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
A   = -5e-2;
Ctr = fem.get('Center0');
dV  = A*invlerp(Ctr(:,1),-15,15)*sin(30*fem.Time);
end