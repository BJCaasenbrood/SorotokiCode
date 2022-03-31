clr;
%% mesh
%sdf = SDF(10,10);
H = 20;
W = 15;

[V,F] = VertexSet(W,H);

msh = Mesh(V,F);
msh = msh.generate();

msh.show();

fem = Fem(msh,'TimeStep',1/1250,'TimeEnd',1,...
    'BdBox',[-5*W,5*W,-15,H],'Linestyle','none');

fem.Material = Dragonskin10(10);

fem = fem.AddConstraint('Gravity',[],[0,-9.81e3]);
fem = fem.AddConstraint('Contact',@(x) SDF(x,W),[0,0]);
id = fem.FindNodes('NE');

%fem = fem.AddConstraint('Spring',id,[0,-1e3]);
%fem = fem.AddConstraint('Pressure',{[9,12]},@(x) -2e-3*sin(30*x.Time));
%fem = fem.AddConstraint('Pressure',{[2,3]},@(x)  -2e-3*sin(40*x.Time));
fem = fem.AddConstraint('Load',fem.FindNodes('Top'),[2e-5,0]);

fem = fem.simulate();

%% movie
t = fem.Log.t; close all;
figure(101);


for ii = 1:fps(t,60):numel(t)
    fem.set('Node',fem.Log.Node{ii});
    fem.show();
    %caxis([0 1e-5]);
    axis([-60 130 -120 30]);
    drawnow();
end


function [V,F] = VertexSet(W,H)

V = [0,0;
     W/3,0;
     W/3,H/3;
     0,H/3;
     W/3,H
     0,H;
     2*W/3,H/3;
     2*W/3,H;
     2*W/3,0;
     W,0;
     W,H/3;
     W,H];
 
V(:,1) = V(:,1) -  W/2; 
 
F = [1,2,3,4;
     4,3,5,6
     3,7,8,5
     9,10,11,7;
     11,7,12,8];

end

function D = SDF(x,W)
    %S = sRectangle(-W,2*W,-150,-5);
    S = sLine(2*W,-W,-2,-5.5);
    D = S.eval(x);
end