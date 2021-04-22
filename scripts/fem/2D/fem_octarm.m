clr;
%% generate mesh
Simp  = 0.001;
GrowH = 1;
MinH  = 0.15;
MaxH  = 2.5;

msh = Mesh('Octarm.png','BdBox',[0,120,0,12],'SimplifyTol',Simp,...
    'Hmesh',[GrowH,MinH,MaxH]);

msh = msh.show();

%% generate fem model
fem = Fem(msh);
fem = fem.set('TimeStep',1/10,'Linestyle','none');

%% add boundary constraint
CP1 = [60,10.5];  % control point 1
CP2 = [60,1.5];   % control point 1
CP3 = [120,8.7];  % control point 1
CP4 = [120,3.3];  % control point 1
F1  = 15e-3;      % control force 1
F2  = 0;          % control force 1
F4  = -5e-3;      % control force 1
F3  = 0;          % control force 1

fem = fem.AddConstraint('Support',fem.FindNodes('Left'),[1,1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Location',CP1),[0,F1]);
fem = fem.AddConstraint('Load',fem.FindNodes('Location',CP2),[0,F2]);
fem = fem.AddConstraint('Load',fem.FindNodes('Location',CP3),[0,F3]);
fem = fem.AddConstraint('Load',fem.FindNodes('Location',CP4),[F4,0]);

%% assign material
fem.Material = Dragonskin10;

%% solve
f = figure(101);
[~,tmpNode] = fem.solve();

b = uicontrol('Parent',f,'Style','slider','Position',[81,24,419,23],...
              'value',1,'min',1, 'max',length(tmpNode));
          
b.Callback = @(s,e) updateNode(s,e,tmpNode,[0 120 -50 50]);

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

