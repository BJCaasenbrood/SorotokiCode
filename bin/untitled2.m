clr;

C1 = sCube([-2,2,-5,5,0,5]);
C2 = sCube([-1,1,-5,4,0,4]);

% D = C1 - C2;
% % %D.show();
% % 
% % obj = Gmodel(D,'Quality',30,'Shading','Face');
% % 
% % obj.bake().render();

D = sTorus(0,0,0,5,2);
%D = sSphere(1);
%S1 = sCube(-1,1,0,1,0,1);

%tic;
%obj = Gmodel(D,'Quality',32);
%toc;

%obj.bake().render();
%tic;
%tic;

[F,V] = DualContour(D,32);

%toc;

patch('Faces',F,'Vertices',V,'FaceColor','w');

% tic;
obj = Gmodel(V,F);
% toc;
% 
% % % % 
obj.bake().render();