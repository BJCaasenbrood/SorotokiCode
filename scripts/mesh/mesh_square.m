clc; clear; close all;
%% set signed distance function
R = 1;

sdf = @(x) dRectangle(x,0,R,0,R);

%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[0,R,0,R],...
              'NElem',150,...
              'MaxIteration',500,...
              'ShowMeshing',true);
      
msh = msh.generateMesh;