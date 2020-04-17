clc; clear; close all;
%% set signed distance function
R = 1;

sdf = @(x) dCircle(x,0,0,R);

%% generate mesh
msh = Mesh(sdf);
msh = msh.set('BdBox',[-R,R,-R,R],...
              'NElem',150,...
              'MaxIteration',500,...
              'ShowMeshing',true);
      
msh = msh.generate();