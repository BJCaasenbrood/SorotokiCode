clr;
%% model
obj = Gmodel('Bunny.stl');

%% set texture
obj.Texture = grey;
obj.Emission = [0.75 0.75 0.75];

% obj.set('SSSPower',1.6,...
%         'AORadius',0.5,...
%         'AOPower',2,...
%         'AO',true);
% 
obj.bake();
    
%% set view
obj.render(); 
obj.ground();
