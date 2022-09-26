function rebuildMexSorotoki()
Path = cdsoro;
cout('Text', '\n');
cout('Text','* Proceeding with MEX libary\n');
Request = input(['* It is highly recommended to use MEX files,',...
    ' do you want to build MEX? (y/n)'],'s');

switch(Request)
    case('y'); buildMexbool = 1;
    case('n'); buildMexbool = 0;
    otherwise; buildMexbool = 0;
end

if buildMexbool
   cout('Text','* Building MEX executables...\n');
   cd([Path,'/src/fem/mex/polardecomposition']);
   cout('Text','* Finding polar decomposition folder...\n');
   cout('Keyword','* building PolarDecompositionFast.mex...\n');
   cout('Green','\t ');
   generateMexPD;
   
   cd([Path,'/src/fem/mex/nonlinearstrainmat']);
   cout('Text','* Finding nonlinear strain mat folder...\n');
   cout('Keyword','* building NonlinearStrainOperatorFast.mex...\n');
   cout('Green','\t ');
   generateMexNSOF;
   
   cd([Path,'/src/fem/mex/localsNH']);
   cout('Text','* Finding locals NeoHookean folder...\n');
   cout('Keyword','* building LocalsNHFast.mex...\n');
   cout('Green','\t ');
   generateMexLNH;
   
   cd([Path,'/src/fem/mex/localsYH']);
   cout('Text','* Finding locals Yeoh folder...\n');
   cout('Keyword','* building LocalsYHFast.mex...\n');
   cout('Green','\t ');
   generateMexLYH;
   
   cd([Path,'/src/fem/mex/localsMN']);
   cout('Text','* Finding locals Mooney folder...\n');
   cout('Keyword','* building LocalsMNFast.mex...\n');
   cout('Green','\t ');
   generateMexLMN;
   
   cd([Path,'/src/model/mex/forwardkin']);
   cout('Text','* Finding forward kinematics folder...\n');
   cout('Keyword','* building computeForwardKinematicsFast.mex...\n');
   cout('Green','\t ');
   generateMexFKF;
   
   cd([Path,'/src/model/mex/lagrangian']);
   cout('Text','* Finding lagrangian computation folder...\n');
   cout('Keyword','* building computeLagrangianFast.mex...\n');
   cout('Green','\t ');
   generateMexCLF;
   

%    cd([Path,'/src/gmodel/mex/marchingcube']);
%    cout('Text','* Finding marching cube folder...\n');
%    cout('Keyword','* building MarchingCubesFast.mex...\n');
%    cout('Green','\t ');
%    generateMexMCF;
%    cout('Text','* MEX build completed!\n');
   
   cd([Path,'/src/gmodel/mex/trinormal']);
   cout('Text','* Finding trinormal folder...\n');
   cout('Keyword','* building TriangleNormalFast.mex...\n');
   cout('Green','\t ');
   generateMexTNF;
   cout('Text','* MEX build completed!\n');
   cdsoro;
else
   cout('Error',['* Not building MEX might result in ,',...
    ' incompatability in future release of SOROTOKI!']);
end
end