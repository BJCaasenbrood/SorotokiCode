clr;
%% model
obj = Gmodel(SDF,'Quality',100);

%% set texture
obj.set('Texture',softmath);
obj.bake().render(); 
   
view(-30,30);
obj.ground();
obj.update();

function y = SDF
    C1 = sCube(0,1,0,1,0,1);
    S1 = sSphere(0,0,1,.5);
    S2 = sSphere(0,0,0.5,1);   
    y = (C1 - S1)/S2;
end