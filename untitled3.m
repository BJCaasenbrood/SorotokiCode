clr;
file1 = 'Cube0.stl';
% file2 = 'Fiber_stub.stl';
% file3 = 'Fiber_stub_chamfer.stl';
% file4 = 'Fiber_stub_fillet.stl';
% 
% 
% obj1 = Gmodel(file2); obj1.set('Shading','Face');
% 
% obj1.render().ground();
model = createpde(1);
importGeometry(model,file1);

generateMesh(model,'Hmax',0.3)

pdeplot3D(model)