clr;
%% generate mesh
Simp  = 0.02;
GrowH = 1.5;
MinH  = 3;
MaxH  = 5;

msh = Mesh('FEMbot.png','BdBox',[0,170,0,100],'SimplifyTol',Simp,...
    'Hmesh',[GrowH,MinH,MaxH]);

msh = msh.show();

%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/40,'Linestyle','none','SolverPlotType','Svm');

%% assign material
fem.Material = NeoHookeanMaterial;

%% add boundary constraint
CP1 = [20,48];  % control point 1
CP2 = [150,48]; % control point 2
F1  = -5;     % control force 1
F2  = 0;      % control force 2

fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[1,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Location',CP1),[F1,0]);
fem = fem.AddConstraint('Load',fem.FindNodes('Location',CP2),[F2,0]);

%% solve
f = figure(101);
[~,tmp1] = fem.solve();

%% reverse-boundary solve
fem = fem.reset();
fem = fem.AddConstraint('Support',fem.FindNodes('Bottom'),[1,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Location',CP1),[-F1,0]);
fem = fem.AddConstraint('Load',fem.FindNodes('Location',CP2),[F2,0]);

[~,tmp2] = fem.solve();

tmpNode = [tmp1(end),tmp2(end)];%[fliplr(tmp1(2:end)),tmp2];

b = uicontrol('Parent',f,'Style','slider','Position',[81,24,419,23],...
              'value',length(tmpNode)/2,'min',1, 'max',length(tmpNode));
          
b.Callback = @(s,e) updateNode(s,e,tmpNode,[-20 180 -20 130]);

function updateNode(src,~,Nd,Axis)
    class = whoClasses('Fem');
    
    low = floor(src.Value);
    upp = ceil(src.Value);

    N = lerp(Nd{low},Nd{upp},src.Value - low);
    
    for i = 1:length(class)
        class{i}.set('Node',N);
        class{i}.show('Un');
    end
    axis(Axis);
end

