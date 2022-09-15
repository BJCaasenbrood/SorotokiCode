clr;
%% mesh

sdf = sCircle(18,25,5) - sCircle(18,25,3.75);
msh = Mesh(sdf,'NElem',50);
msh = msh.generate();

%%
fem = Fem(msh,'TimeStep',1/1e3,'TimeEnd',1.5,'BdBox',[-35,35,-15,30],...
    'Linestyle','none');

fem.Material = NeoHookeanMaterial(0.02,0.4);
fem.Material.Cfr = 3e-6;

fem = fem.addGravity();
fem = fem.addContact(SDF(15),[0,0]);

fem = fem.simulate();

%% movie
t = fem.Log.t; close all;
figure(101); clf;

for ii = 1:fps(t,200):numel(t)
    fem.set('Node',fem.Log.Node{ii});
    
    fem.show();
    caxis([-1e-6 1e-6])
    drawnow();
    background(metropolis);
    
        
    if ii == 1
        gif('soft_bounce.gif','frame',gcf,'nodither');
    else
        gif; 
    end

end

%% Signed distance environment
function D = SDF(W)
    S1 = sLine(4,-W,-10,-5);
    S2 = sLine(W,4,0,-10);
    S3 = sLine(-W,-W,-1,1);
    D = S1 + S2 + S3;
end