
%% mesh
sdf = sCircle(18,25,5) - sCircle(18,25,3.5);
msh = Mesh(sdf,'NElem',55);
msh = msh.generate();

%% fem model
fem = Fem(msh,'TimeStep',1/1200);

fem.Material = NeoHookeanMaterial(0.01,0.49);
fem.Material.Rho = fem.Material.Rho;
fem.Material.Cfr  = 0.5;
fem.Material.Rr   = 5;

fem = fem.addGravity();
fem = fem.addContact(SDF(15),[0,0]);

fig(101,[9.5,9.5]);
fem = fem.simulate();

%% Signed distance environment
function D = SDF(W)
    S1 = sLine(4,-W,-10,-5);
    S2 = sLine(W,4,0,-10);
    S3 = sLine(-W,-W,-1,1);
    D = S1 + S2 + S3;
end