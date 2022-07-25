clr;
%% mesh

sdf = sCircle(18,25,5) - sCircle(18,25,3.75);
msh = Mesh(sdf,'NElem',60);
msh = msh.generate();

%%
fem = Fem(msh,'TimeStep',1/800,'TimeEnd',1.5,'BdBox',[-35,35,-15,30],...
    'Linestyle','-');

fem.Material = NeoHookeanMaterial(0.1,0.4);
fem.Material.Rho  = 5*fem.Material.Rho;
fem.Material.Cfr  = 1e-4;

fem = fem.AddConstraint('Gravity',[],[0,-9.81e3]);
fem = fem.AddConstraint('Contact',@(x) SDF(x,15),[0,0]);

fem = fem.simulate();

%% movie
t = fem.Log.t; close all;
figure(101); clf;

for ii = 1:fps(t,200):numel(t)
    fem.set('Node',fem.Log.Node{ii});
    
    fem.show();
    caxis([-1e-6 1e-6])
    drawnow();
    background(gitpage);
%     if ii == 1, gif('fem_bouncingball.gif',...
%             'frame',gcf,'nodither','DelayTime',1/30);
%     else, gif; 
%     end
end

%% Signed distance environment
function D = SDF(x,W)
    S1 = sLine(4,-W,-10,-5);
    S2 = sLine(W,4,0,-10);
    S3 = sLine(-W,-W,-1,1);
    S = S1 + S2 + S3;
    D = S.eval(x);
end