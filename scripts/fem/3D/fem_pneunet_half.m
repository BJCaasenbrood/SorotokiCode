clr;
%% meshing
msh = Mesh('PneuNetHalf.stl','Hmesh',[2,3,5]);
%msh = msh.generate();

%% fem
% fem = Fem(msh,'TimeStep',1/30,'Linestyle','-','TimeEnd',3);
% fem.Material = Ecoflex0030(30);
% 
% %%
% fem = fem.addSupport(fem.FindNodes('Left'),[1,1,1]);
% fem = fem.addSupport(fem.FindNodes('Front'),[0,1,0]);
% fem = fem.addGravity([0,0,-9810]);
% 
% %%
% fem.solve();
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