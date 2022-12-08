clr;
%% meshing
msh = Mesh('PneuNetHalf.stl','Hmesh',[1,5,6]);
msh = msh.generate();

%% fem
P = 5*kpa;

fem = Fem(msh,'TimeStep',1/250,'TimeEnd',3,'ResidualNorm',Inf);
fem.Material = Dragonskin10();%NeoHookeanMaterial(0.5,0.3);%Ecoflex0030(30);

fem = fem.addSupport(fem.FindNodes('Left'),[1,1,1]);
fem = fem.addSupport(fem.FindNodes('Front'),[0,1,0]);

%fem = fem.addGravity();

% fem = fem.addTendon(fem.FindNodes('Bottom'),[0,0,-1e-3]);
% fem = fem.addTendon(fem.FindNodes('Right'),[1e-3,0,0]);

%%
%msh.show();

 pln = [5.14,5.14,-12,-2,3.75,18.75];
 fem = fem.addTendon(fem.FindNodes('Box',pln),[-150*P,0,0]);
 pln = [9,9,-12,-2,3.75,18.75];
 fem = fem.addTendon(fem.FindNodes('Box',pln),[150*P,0,0]);
 pln = [108.3,108.3,-12,-2,3.75,18.75];
 fem = fem.addTendon(fem.FindNodes('Box',pln),[-150*P,0,0]);
 pln = [111.55,111.55,-12,-2,3.75,18.75];
 fem = fem.addTendon(fem.FindNodes('Box',pln),[150*P,0,0]);
 pln = [108.32,111.55,-12,-2,3.75,3.75];
 fem = fem.addTendon(fem.FindNodes('Box',pln),[0,0,30*P]);
 pln = [108.32,111.55,-12,-2,18.75,18.75];
 fem = fem.addTendon(fem.FindNodes('Box',pln),[0,0,-30*P]);

% 
for ii = 1:10
    DX = 10.2;
    pln = [16.0748+(ii-1)*DX,16.0748+(ii-1)*DX,...
        -12,-2,3.75,18.75];
    fem = fem.addTendon(fem.FindNodes('Box',pln),[-150*P,0,0]);
end

for ii = 1:10
    DX = 10.2;
    pln = [19.45+(ii-1)*DX,19.45+(ii-1)*DX,...
        -12,-2,3.75,18.75];
    fem = fem.addTendon(fem.FindNodes('Box',pln),[150*P,0,0]);
end

for ii = 1:10
    DX = 10.2;
    pln = [16.0748+(ii-1)*DX,19.45+(ii-1)*DX,...
        -12,-2,3.75,3.75];
    fem = fem.addTendon(fem.FindNodes('Box',pln),[0,0,30*P]);
end

for ii = 1:10
    DX = 10.2;
    pln = [16.0748+(ii-1)*DX,19.45+(ii-1)*DX,...
        -12,-2,18.75,18.75];
    fem = fem.addTendon(fem.FindNodes('Box',pln),[0,0,-30*P]);
end

% %fem = fem.addSupport(fem.FindNodes('Front'),[0,1,0]);
% id  = fem.FindNodes('Box',[60.4,60.4,-12,-2,3.75,18.75]);
% fem = fem.addTendon(id,[P,0,0]);
% %fem = fem.addLoad(fem.FindNodes('Right'),[0,0,-1e-4]);
% 
% id = fem.FindNodes('Box',[56.7,56.7,-12,-2,3.75,18.75]);
% fem = fem.addTendon(id,[-P,0,0]);
%fem = fem.addGravity([0,0,-1]);
% 
% 
% id = fem.FindNodes('Bottom');
% fem = fem.addTendon(id,[0,0,-P]);

%%
fem.solve();
% 
% %%
% close all;
% figure(101);
% 
% t = fem.Log.t;
% 
% for ii = 1:fps(t,20):numel(t)
%     fem.set('Node',fem.Log.Node{ii});
%     fem.show();
%     axis([-80 120 -30 30 -130 50]);
%         
%     drawnow();
% end