clr;
%% mesh
Y0 = 25;

sdf = sCircle(18,Y0,5) - sCircle(18,Y0,3.5);
msh = Mesh(sdf,'NElem',50);
msh = msh.generate();

%% f%% fem model
fem = Fem(msh,'TimeStep',1/800,'TimeEnd',1.5,...
    'BdBox',[-35,35,-15,30],'Linestyle','-');

fem.Material = NeoHookeanMaterial(0.1,0.4);
fem.Material.Rho = fem.Material.Rho;
fem.Material.Zeta = 1e-6;
fem.Material.Cfr  = 0.5;
fem.Material.Rr   = 5;

fem = fem.addGravity();
fem = fem.addContact(SDF(15),[0,0]);

fem = fem.simulate();

%% Signed distance environment
function D = SDF(W)
    S1 = sLine(4,-W,-10,-5);
    S2 = sLine(W,4,0,-10);
    S3 = sLine(-W,-W,-1,1);
    D = S1 + S2 + S3;
    %D = sLine(0,-1,0,0);
end