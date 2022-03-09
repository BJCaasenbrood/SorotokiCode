N    = rand(4,1);
dNdx = rand(4,2);
F = eye(3);

[Bn,Bg,NN,tau] = NonlinearStrainOperatorFast(N,dNdx,F)