clr;
%% loading .stl file
obj = Gmodel('Bunny.stl');

%% set material texture
obj.Texture = sunshine;
obj.Emission = [0.85 0.85 0.9];
%  
% obj.set('SSSPower',0.75,'SSSRadius',0.3,'SSS',true,...
%         'AOPower',1.85,'AORadius',0.3,'AO',true);
%      
%% rendering    
obj.bake().render(); 
obj.ground();
