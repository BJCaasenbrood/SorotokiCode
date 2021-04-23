clr;
%% generate signed distance functions (i.e., circles)
f1 = @(x) sqrt(x(:,1).^2 + x(:,2).^2)- 0.5;
f2 = @(x) sqrt((x(:,1)+0.5).^2 + x(:,2).^2) - 0.5;

%% generate SDF class
sdf1 = Sdf(f1);
sdf2 = Sdf(f2);

% assigning bounding box
sdf1.BdBox = [-1,1,-1,1];
sdf2.BdBox = [-1.5,0.5,-1,1];

%% math operations
f = sdf1-sdf2;
f.show(); view(0,90);
pause(1);

%% generate mesh from SDF
msh = Mesh(f,'NElem',1e2);
msh = msh.generate();
msh.show();