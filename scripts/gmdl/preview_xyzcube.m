clr;
%% loading .stl file
% obj = Gmodel('Cube0.stl');
% 
% %% set material texture
% obj.Texture = prusa;
% obj.Emission = 1.05*[0.95 0.05 0.05];
%  
% obj.set('SSSPower',0.4,'SSSRadius',0.3,'SSS',true,...
%         'AOPower',2.85,'AORadius',0.3,'AO',true);
%      
% %% rendering    
% obj.bake().render(); 
% obj.ground();

[V,F,TF,VT] = dimcube(0.1);

cla; patcht(F,V,TF,VT,xyzcube);
%patch('faces',TF,'vertices',VT,'facecolor','none');