clr; 

q = [12,0,0]*kpa;

Y = pccspace(100,1); 
[base,rig,grip] = rig_softmanipulator(Y); 

G = rig.FK(q);

rig     = rig.computeFK(q); 
grip.g0 = G(:,:,end);
grip = grip.computeFK(-0.03);


grip.update;
rig.update;

axis tight;