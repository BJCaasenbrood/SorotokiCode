clc; clear; close all;
%% set signed distance function
sdf = @(x) dRectangle(x,0,1,0,1);

%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[0,1,0,1],...
              'NElem',50,...
              'ShowMeshing',true);
      
msh = msh.generateMesh;