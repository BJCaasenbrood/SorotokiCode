clr;
% generate mesh from sdf
sdf = Sdf( @(x) TensileBone(x,11,2,3,1,1), 'BdBox', [0,10,0,10]);

msh = Mesh(sdf,'NElem',150);
msh = msh.generate();

% generate fem model from mesh
fem = Fem(msh,'TimeStep',1/75,'LineStyle','none');

% add boundary conditions
fem = fem.addSupport('sw',[1,1]);
fem = fem.addSupport('bottom',[1,1]);
fem = fem.addDisplace('top',[0,25]);

% assign material
fem = fem.addMaterial( NeoHookean(0.1,0.3) );
% fem = fem.addMaterial( Yeoh([0.1, 0.01, 0],[100,100,100]) );

% solving
fem.solve();

%% Tensile bone SDF
function D = TensileBone(P,H,W,T,D,R)
    dD = 0.5*(W-D);
    dT = 0.5*(H-T);

    R1 = dRectangle(P,0,W,0,H);
    R2 = dRectangle(P,0,dD,dT,dT+T);
    R3 = dRectangle(P,W-dD,W,dT,dT+T);
    C1 = dCircle(P,dD-R,dT,R);
    C2 = dCircle(P,dD-R,dT+T,R);
    C3 = dCircle(P,W-dD+R,dT,R);
    C4 = dCircle(P,W-dD+R,dT+T,R);
    D0 = dDiff(dDiff(dDiff(R1,R2),C1),C2);
    D = dDiff(dDiff(dDiff(D0,R3),C3),C4);
end