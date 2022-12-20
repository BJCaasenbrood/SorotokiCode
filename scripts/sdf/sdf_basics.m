clr;
%% render
fig(101,[9.5,9.5])
S = sCylinder(2);
obj = Gmodel(S,'Quality',120);

view(30,40);
obj.bake.render();

export_fig render.png

