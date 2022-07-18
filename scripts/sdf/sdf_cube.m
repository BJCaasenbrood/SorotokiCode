clr;
%%
sdf0 = sCube(1);

%% render equivalent 
sdf = sRectangle(1);
sdf.BdBox = [-2,2,-2,2];

sdf.show(); hold on;
sdf0.skeleton();