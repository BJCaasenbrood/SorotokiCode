clr;
H1 = Gmodel('FrameSRMHolder.stl','Shading','Face'); 
H1 = Blender(H1,'Scale',1/64);
L1 = Gmodel('Pneulink.stl');    
C1 = Gmodel('ConnectorRedux.stl','Shading','Face');
C2 = Gmodel('ConnectorRedux.stl','Shading','Face');
C2 = Blender(C2,'Scale',{'z',0.4});
G1 = Gmodel('GripperFingerRedux.stl');

H1.Texture = egg;
L1.Texture = bluebase;
G1.Texture = L1.Texture;
C1.Texture = H1.Texture;

[L1, He] = Blender(L1,'Curve',{'PCC+',90,-90,1.5});

G1 = Blender(G1,'Curve',{'PCC',0,-39});
G1 = Blender(G1,'Rotate',{'x',-30});

C2 = Blender(C2,'Translate',{'z',0.1});
G1 = Blender(G1,'Translate',{'z',0.15});
G1 = Blender(G1,'Translate',{'y',0.25});
G1 = Blender(G1,'CArray',{'z',3});

G1 = Blender(G1,'SE3',He);
C1 = Blender(C1,'SE3',He);
C2 = Blender(C2,'SE3',He);

G1 = Blender(G1,'Rotate',{'y',180});
L1 = Blender(L1,'Rotate',{'y',180});
C1 = Blender(C1,'Rotate',{'y',180});
C2 = Blender(C2,'Rotate',{'y',180});
H1 = Blender(H1,'Rotate',{'y',180});

L1.bake.render();
C1.bake.render();
C2.bake.render();
G1.bake.render();
H1.bake.render();

cyl = Gmodel('Cylinder.stl'); cyl.Texture = diffuse(0.03)*1.1;
cyl = Blender(cyl,'Scale',{'z',0.35});
cyl = Blender(cyl,'Scale',{'xy',0.65});
cyl = Blender(cyl,'Translate',{'z',0.4});
cyl = Blender(cyl,'SE3',He);
cyl = Blender(cyl,'Rotate',{'y',180});

cyl.render(); 
view(50,30);



axis tight;
view(190,10);