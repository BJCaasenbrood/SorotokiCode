function [Bn,Bg,NN,tau] = NonlinearStrainOperatorFast(N,dNdx,F)
nn = length(N);
mm = size(dNdx,2);
zz = mm*nn;

NN = zeros(mm,zz);

id1 = 1:mm:zz;
id2 = 2:mm:zz;

NN(1,id1) = N.';
NN(2,id2) = N.';
dNdxX = dNdx(:,1).';
dNdxY = dNdx(:,2).';

if mm == 3
    NN(3,3:mm:zz) = N.';
    dNdxZ = dNdx(:,3).';
else
    dNdxZ = 0;
end

Bn = zeros((mm-1)*3,zz);
Bg = zeros((mm-1)*4+(mm-2),zz);

if mm == 2 % 2-dimensional
    Bn(1,id1) = dNdxX*F(1,1);
    Bn(1,id2) = dNdxX*F(2,1);
    Bn(2,id1) = dNdxY*F(1,2);
    Bn(2,id2) = dNdxY*F(2,2);
    Bn(3,id1) = dNdxX*F(1,2) + dNdxY*F(1,1);
    Bn(3,id2) = dNdxX*F(2,2) + dNdxY*F(2,1);
    
    Bg(1,id1) = dNdxX;
    Bg(2,id1) = dNdxY;
    Bg(3,id2) = dNdxX;
    Bg(4,id2) = dNdxY;
    
else % 3-dimensional
    Bn(1,1:mm:mm*nn) = dNdxX*F(1,1);
    Bn(1,2:mm:mm*nn) = dNdxX*F(2,1);
    Bn(1,3:mm:mm*nn) = dNdxX*F(3,1);
    Bn(2,1:mm:mm*nn) = dNdxY*F(1,2);
    Bn(2,2:mm:mm*nn) = dNdxY*F(2,2);
    Bn(2,3:mm:mm*nn) = dNdxY*F(3,2);
    Bn(3,1:mm:mm*nn) = dNdxZ*F(1,3);
    Bn(3,2:mm:mm*nn) = dNdxZ*F(2,3);
    Bn(3,3:mm:mm*nn) = dNdxZ*F(3,3);
    Bn(4,1:mm:mm*nn) = dNdxX*F(1,2) + dNdxY*F(1,1);
    Bn(4,2:mm:mm*nn) = dNdxX*F(2,2) + dNdxY*F(2,1);
    Bn(4,3:mm:mm*nn) = dNdxX*F(3,2) + dNdxY*F(3,1);
    Bn(5,1:mm:mm*nn) = dNdxY*F(1,3) + dNdxZ*F(1,2);
    Bn(5,2:mm:mm*nn) = dNdxY*F(2,3) + dNdxZ*F(2,2);
    Bn(5,3:mm:mm*nn) = dNdxY*F(3,3) + dNdxZ*F(3,2);
    Bn(6,1:mm:mm*nn) = dNdxX*F(1,3) + dNdxZ*F(1,1);
    Bn(6,2:mm:mm*nn) = dNdxX*F(2,3) + dNdxZ*F(2,1);
    Bn(6,3:mm:mm*nn) = dNdxX*F(3,3) + dNdxZ*F(3,1);
    
    Bg(1,1:mm:mm*nn) = dNdxX;
    Bg(2,1:mm:mm*nn) = dNdxY;
    Bg(3,1:mm:mm*nn) = dNdxZ;
    Bg(4,2:mm:mm*nn) = dNdxX;
    Bg(5,2:mm:mm*nn) = dNdxY;    
    Bg(6,2:mm:mm*nn) = dNdxZ;   
    Bg(7,3:mm:mm*nn) = dNdxX;
    Bg(8,3:mm:mm*nn) = dNdxY;    
    Bg(9,3:mm:mm*nn) = dNdxZ;     
end

tau = 1;

end 