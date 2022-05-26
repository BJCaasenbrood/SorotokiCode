clr;
%%
sdf0 = sCube(1);
% obj = Gmodel(,'Quality',60,...
%     'Texture',grey,'ShowProcess',0);
% obj.bake.render();

%% render equivalent 
sdf = sRectangle(1);
sdf.BdBox = [-2,2,-2,2];

sdf.show(); hold on;
sdf0.skeleton();