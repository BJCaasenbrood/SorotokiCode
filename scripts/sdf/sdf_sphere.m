clr;
%%
S = sSphere(1)
obj = Gmodel(S,...
    'Texture',grey,'ShowProcess',0);
obj.bake.render();

%% render equivalent 
sdf = sCircle(1);
sdf.BdBox = [-2,2,-2,2];

hold on;
sdf.show();
S.skeleton();