% clr;
% %% generate SDF class
% R1 = sCube(0,100,-10,10,0,4);
% R2 = sCube(5,20,-10,10,0,20);
% %sdf2 = sCircle(0.5, 0, 0.5);
% 
% %% math operations
% S = R1;
% W = 13 + 2;
% 
% for ii = 0:5
%     S = S + sCube(5+W*ii,18+W*ii,-10,10,0,20);
% end
% 
% obj = Gmodel(S,'Quality',45);
% obj.bake.render();
