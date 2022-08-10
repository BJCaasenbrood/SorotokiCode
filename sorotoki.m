function sorotoki(arg)
%SOROTOKI is the installer function of SOROTOKI.
%   SOROTOKI() can be executed in the command window and it will start the
%   installation aumatically. Make sure you have an active internet
%   connection to start the installation, which will be compared to the
%   most current (stable) version of the toolkit. Some embedded functions
%   for the SOROTOKI installer are:
%
%   sorotoki install       % installer
%   sorotoki -i            % ...
%   sorotoki check         % check installation/toolkit
%   sorotoki -c            % check installation/toolkit
%   sorotoki demo          % lists demos
%   sorotoki -d            % lists demos
%   sorotoki help          % help
%   sorotoki -h            % help
%   sorotoki version       % version
%   sorotoki -v            % version
%
%   Also, some additional toolboxes are required:
%
%   - Image Processing Toolbox
%   - Partial Differential Equation Toolbox
%   - MATLAB Coder
%   - Image Processing Toolbox
%
%   The installer will automaticall check for these. Furthermore, it will
%   also prompt for configurating a startup file. This will ensure that
%   sorotoki is also loaded when opening MATLAB, and can be found under
%   .../MATLAB/startup.m in the MATLAB's main folder on your desktop.
%
%   Last edit: Jul 15, 2022.

if nargin < 1, clc; clear; arg = 'install'; end
% ------------------------------------------------------------------------
switch(arg)
    case('install');   setupToolkit;
    case('uninstall'); unloadSorotoki;
    case('remove');    unloadSorotoki;
    case('demo');      showDemo;

    case('cd');        getPath;
    case('check');     verifySorotoki;
    case('version');   disp(soropatch(1));
    case('--version'); disp(soropatch(1));
    case('-v');        disp(soropatch(1));
    case('-c');        verifySorotoki;
    case('-i');        setupToolkit;
    case('-h');        help sorotoki;
    case('-d');        showDemo;
    case('-u');        unloadSorotoki;
end

end

% -------------------------------------------------------------- show demos
function showDemo
clr;
set = {'Meshing a 2D-circle',...
       'Rendering stanford bunny',...
       'Nonlinear finite element',...
       '3D nonlinear finite element',...
       'Topology optimization of Pneu-net',...
       'Modeling a Cosserat beam',...
       'Closed-loop control soft arm'};

request = action(set);

switch(request)
    case(1); open mesh_circle;
    case(2); open preview_bunny;
    case(3); open fem_bellow;
    case(4); open fem_twisting_beam;
    case(5); open opt_pneunet;
    case(6); open mdl_free_motion;
    case(7); open mdl_closedloop_setpoint
    otherwise; warning('Please select a demo from the list above.');
end

end
% ----------------------------------------------------------- setup toolkit
function setupToolkit  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all; warning off; beep off;

addpath('src');
addpath('src/__base__');
addpath('src/__base__/auxiliary');
addpath('src/__base__/thirdparty');
addpath('src/__base__/colormap');
addpath('src/__base__/color');
addpath('src/__base__/userinput');
addpath('src/__base__/plot');

skipUpdate = false;
Path = cd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
displayLogoSorotoki;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cout(['* Thank you for installing SOROTOKI, we will',...
    ' start the installation shortly \n']); pause(1.0);

cout('* Starting installation - SOROTOKI');
verFolder = 'src/__version__';
if ~exist(verFolder, 'dir')
mkdir(verFolder);
cout(['* Created directory ',verFolder, ' \n']);
cout(['* Directory ',verFolder, ' added to path \n']);
else
cout(['* Directory ',verFolder, ' added to path \n']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addPath([Path,verFolder]);
global FIT FID
InstallerFile = [Path,'/src/__version__/install.log'];

if exist(InstallerFile,'file')
    delete(InstallerFile); 
end

FIT = fopen(InstallerFile,'w');
fprintf(FIT,'%% started sorotoki installer on ');
fprintf(FIT,[datestr(datetime('today')),'\n']);
fprintf(FIT,'%% install directory = ');
fprintf(FIT,[Path,'\n']);

if ~pingserver
    cout('err','No internet connection! '); 
    str = action({'(y)es, continue without connection',...
            '(n)o, stop installation'},'s');
    if strcmp(str,'N'), return; end
    skipUpdate = true;
end

if ~skipUpdate
cout('* Getting version_file_check from Git repository src/soropatch.m \n');
url = ['https://raw.githubusercontent.com/BJCaasenbrood/',...
    'SorotokiCode/master/src/'];
filename = [verFolder,'/soropatch.m'];
%websave(filename,url);
cout(['* Succesfully downloaded latest patchnotes -- ', filename, '\n']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
libs(1) = includeBase(Path,0);
libs(2) = includeGraphicsModel(Path,0);
libs(3) = includeMesh(Path,0);
libs(4) = includeFiniteElement(Path,0);
libs(5) = includeDynamicModel(Path,0);
libs(6) = includeControl(Path,0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if min(libs) == 1
cout('\n* Libary check completed - all libaries are up-to-date - \n');
else
cout('err','\n* Libary check completed - some libaries are outdated!\n');
cout('* Get the latest version at: ');

fprintf(['<a href="https://github.com/BJCaasenbrood/SorotokiCode">',...
    'https://github.com/BJCaasenbrood/SorotokiCode</a> \n'])

Request = input('* Do you want to continue the installation? (y/n) ','s');

switch(Request)
    case('y'); cout('green','* Proceeding toolkit installation... \n'); pause(.1);
    case('n'); cout('error','* Installation aborted... \n');
    otherwise; error('Key not recognized...');
end

end
Request = input('* Do you want to create a startup file? (y/n)','s');
bool = 1;
switch(Request)
    case('y'); fprintf('* Proceeding generation startup.m file... \n'); pause(.1);
    case('n'); bool = 0;
    otherwise; bool = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if bool == 1
    %global FID;
    StartUpFile   = [userpath,'/startup.m'];
    delete(StartUpFile);
    
    FID = fopen(StartUpFile,'w');
    
    fprintf('* Assigning fixed search paths to MATLAB: \n');
    
    if libs(1), includeBase(Path,1); end
    if libs(2), includeGraphicsModel(Path,1); end
    if libs(3), includeMesh(Path,1); end
    if libs(4), includeFiniteElement(Path,1); end
    if libs(5), includeDynamicModel(Path,1); end
    if libs(6), includeControl(Path,1); end
    
    Request = input(['\n* Do you want to setup the SorotokiCode folder',...
        ' as your main directory (y/n)'],'s');

    switch(Request)
        case('y'), fprintf(FID,'cdsoro; \n');
    end
    
    fclose('all');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = ver;
ToolboxesInstalled = {S.Name}';
Requisites = {'Image Processing Toolbox';'Partial Differential Equation Toolbox';...
    'MATLAB Coder'; 'Image Processing Toolbox'};
HTMLLink = {'<a href="https://nl.mathworks.com/products/image.html">click here to install</a> \n',...
            '<a href="https://nl.mathworks.com/products/pde.html">click here to install</a> \n'};
memb = ismember(Requisites,ToolboxesInstalled);

cout('Text','* Checking for prerequisit toolboxes...\n'); 
pause(1);
reqCheckList = find(memb == false);

if ~isempty(reqCheckList)
cout('Error', '* Missing some required toolboxes! \n');    

    for ii = 1:length(reqCheckList)
        cout('Error', ['\t ->',Requisites{reqCheckList(ii)},'  '])
        cout('Error',HTMLLink{reqCheckList(ii)});    
    end

    Request = input(['\n* Do you want to continue without the',...
        ' recommended prerequisits? (y/n)'],'s');
    
    switch(Request)
        case('y')
        case('n'), cout('Error','* INSTALLATION TERMINATED! \n'); return;
    end
else
    cout('green','* All prerequisit toolboxes present! \n');   
    pause(0.75);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%check if cmake present
cout('Text','* Checking for c++ compilers...\n'); pause(1);

origLD = getenv('LD_LIBRARY_PATH');
setenv('LD_LIBRARY_PATH','');

[status1,~] = system("gcc --version");
[status2,~] = system("g++ --version");
[status3,~] = system("clang --version");

flag = (status1 ~= 0 && status2 ~= 0 && status3 ~= 0);

if flag
    cout('Error', '* No C++ compiler found! Please make sure you have a C/C++ compiler!\n')
    cout('Error', '* Recommended compiler => gcc:')
    cout('Error', '\t <a href="http://mingw-w64.org/doku.php/start">click here to install</a> \n')
    
    cout('Error',['\n* No C++ compiler implies that no changes can be made to',...
    ' Model.lib (Cosserat Beam model) \n']);
    cout('Error', '* Press ENTER to continue installation\n');
    pause();
else
    cout('green','* C++ compiler founded - ')
    if (status1 == 0), cout('green','gcc '); end
    if (status2 == 0), cout('green','g++ '); end
    if (status3 == 0), cout('green','clang'); end
    cout('\n');
end
pause(0.75);

% cout('Text','* Checking for Cmake...\n'); pause(1);
% 
% [status,~] = system("cmake --version");
% if  status~= 0
%     cout('Error', '* Missing Cmake! \n');    
% else
%     cout('green','* Cmake compiler founded\n')
% end
% pause(0.75);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cout('Text', '\n');
cout('green','* INSTALLATION DONE! \n');
cout('Text', '\n');
cout('* SOROTOKI toolkit is succesfully installed and ready-to-use! \n');

Request = input(['* Do you want to run a check-up if SORTOKI',...
    ' is installed correctly? (y/n)'],'s');

bool = 1;
switch(Request)
    case('y'); fprintf('* Proceeding verify_sorotoki.m file... \n'); pause(.1);
    case('n'); bool = 0;
    otherwise; bool = 0;
end

if bool
   arg = verifySorotoki(); 
   cout('Text', '* Finalizing the verification check-up \n');
   cout('Text', '\n');
   if isempty(vertcat(arg{:,2}))
    cout('green','* VERIFICATION PASSED! \n');
   else
    out = cellfun(@(x) ~isempty(x),arg(:,2));
    cout('error','* VERIFICATION FAILED! \n');
    cout('error','* Error(s) in :');
    errors = arg(out,1);
    for ii = 1:length(errors)
        cout('error',['\t',errors{ii}]); 
        cout('\n');
    end
   end
   cout('Text', '\n');
end

pause(.01);
cout(['* The documentation can be found in doc/SorotokiManual.pdf',...
    '. For more in-\nformation on the Soft Robotics Toolkit, visit',...
    ' the GitHub repository at: \n\n']) 
cout(['<a href="https://bjcaasenbrood.github.io/SorotokiCode">',...
    'https://bjcaasenbrood.github.io/SorotokiCode</a> \n'])
cout('Text', '\n');
cout('Text', '* To get started with demonstrations, type ');
cout('String', 'sorotoki(''demo'') \n');
warning on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
% ---------------------------------------------------------------- get path
function Path = getPath(print)
if nargin < 1; print = false;end
Path = cd;
if print, disp(['* cd: ',Path]); end
end
% ------------------------------------------------------------------ BASICS
function x = includeBase(Path,Request)
global FID FIT
cout(['* Adding SOROTOKI libraries to path, ',...
    'this might take a minute...\n\n']);

if Request == 1
fprintf(FID,['%%!INSTALDIR:',strrep(Path,'\','/'),' \n']);
fprintf(FID,'%% base.lib \n');
write2file([Path,'src']);
write2file([Path,'\src\__base__']);
write2file([Path,'\src\__version__']);
write2file([Path,'\src\__base__\auxiliary']);
write2file([Path,'\src\__base__\thirdparty']);
write2file([Path,'\src\__base__\color']);
write2file([Path,'\src\__base__\color\nasa']);
write2file([Path,'\src\__base__\color\specific']);
write2file([Path,'\src\__base__\color']);
write2file([Path,'\src\__base__\colormap']);
write2file([Path,'\src\__base__\userinput']);
write2file([Path,'\src\__base__\plot']);
write2file([Path,'\src\__base__\units']);
write2file([Path,'\data\']);
write2file([Path,'\data\color']);
write2file([Path,'\data\colormap']);
write2file([Path,'\data\matcap']);
write2file([Path,'\data\matcap\img']);
write2file([Path,'\data\stl']);
write2file([Path,'\data\contours']);
write2file([Path,'\scripts\']);
else
addPath(Path);
addPath([Path,'src']);
addPath([Path,'\src\__version__']);
addPath([Path,'\src\__base__']);
addPath([Path,'\src\__base__\auxiliary']);
addPath([Path,'\src\__base__\thirdparty']);
addPath([Path,'\src\__base__\color']);
addPath([Path,'\src\__base__\color\nasa']);
addPath([Path,'\src\__base__\color\specific']);
addPath([Path,'\src\__base__\colormap']);
addPath([Path,'\src\__base__\userinput']);
addPath([Path,'\src\__base__\plot']);
addPath([Path,'\src\__base__\units']);
addPath([Path,'\data\']);
addPath([Path,'\data\color']);
addPath([Path,'\data\colormap']);
addPath([Path,'\data\matcap']);
addPath([Path,'\data\matcap\img']);
addPath([Path,'\data\stl']);
addPath([Path,'\data\contours']);
addPath([Path,'\scripts\']);
pause(.3);

x = checkLibary('base.lib',@(x) basePathConfirm);

if x
    fprintf(FIT,'%% base.lib installed succesfully! \n');
else
    fprintf(FIT,'%% base.lib missing! \n');
end

end
end
% ---------------------------------------------------------------- GRAPHICS
function x = includeGraphicsModel(Path,Request)
global FID FIT

if Request == 1
fprintf(FID,'%% gmodel.lib \n');
write2file([Path,'\scripts\gmdl\']);
write2file([Path,'\livescripts\gmdl\']);
write2file([Path,'\src\gmodel']);
write2file([Path,'\src\gmodel\tools\']);
write2file([Path,'\src\gmodel\matcap\']);
write2file([Path,'\src\gmodel\mex\']);
write2file([Path,'\src\gmodel\mex\marchingcube']);
write2file([Path,'\src\gmodel\mex\trinormal']);
write2file([Path,'\src\gmodel\matcap\img']);
write2file([Path,'\src\gmodel\tools\thirdparty']);
else
addPath([Path,'\livescripts\gmdl\']);
addPath([Path,'\scripts\gmdl\']);
addPath([Path,'\src\gmodel']);
addPath([Path,'\src\gmodel\tools\']);
addPath([Path,'\src\gmodel\tools\thirdparty']);
addPath([Path,'\src\gmodel\matcap\']);
addPath([Path,'\src\gmodel\mex\']);
addPath([Path,'\src\gmodel\mex\marchingcube']);
addPath([Path,'\src\gmodel\mex\trinormal']);
addPath([Path,'\src\gmodel\matcap\img\']);
addPath([Path,'\src\gmodel\matcap\tools\']);
pause(.3);

x = checkLibary('gmodel.lib',@(x) graphicsmodelPathConfirm);

if x, fprintf(FIT,'%% gmodel.lib installed succesfully! \n');
else, fprintf(FIT,'%% gmodel.lib missing! \n');
end

end
end
% -------------------------------------------------------------------- MESH
function x = includeMesh(Path,Request)
global FID FIT

if Request == 1
fprintf(FID,'%% mesh.lib \n');
write2file([Path,'\src\mesh']);
write2file([Path,'\scripts\mesh\']);
write2file([Path,'\scripts\mesh\SDF']);
write2file([Path,'\src\mesh\tools\']);
write2file([Path,'\src\mesh\shapes\']);
write2file([Path,'\src\mesh\operators\']);
else
addPath([Path,'\scripts\mesh\']);
addPath([Path,'\scripts\mesh\SDF']);
addPath([Path,'\src\mesh']);
addPath([Path,'\src\mesh\tools\']);
addPath([Path,'\src\mesh\shapes\']);
addPath([Path,'\src\mesh\operators\']);
pause(.3);

x = checkLibary('mesh.lib',@(x) meshPathConfirm);

if x
    fprintf(FIT,'%% mesh.lib installed succesfully! \n');
else
    fprintf(FIT,'%% mesh.lib missing! \n');
end

end
end
% --------------------------------------------------------------------- FEM
function x = includeFiniteElement(Path,Request)
global FID FIT

if Request == 1
fprintf(FID,'%% fem.lib \n');
write2file([Path,'\scripts\fem\']);
write2file([Path,'\scripts\fem\2D']);
write2file([Path,'\scripts\fem\3D']);
write2file([Path,'\scripts\opt']);
write2file([Path,'\src\fem']);
write2file([Path,'\src\fem\tools\']);
write2file([Path,'\src\fem\tools\tpswarp']);
write2file([Path,'\src\fem\tools\mma']);
write2file([Path,'\src\fem\mex\']);
write2file([Path,'\src\fem\materials\']);
write2file([Path,'\src\fem\materials\samples']);
write2file([Path,'\src\fem\mex\nonlinearstrainmat']);
write2file([Path,'\src\fem\mex\polardecomposition']);
write2file([Path,'\src\fem\mex\localsYH']);
write2file([Path,'\src\fem\mex\localsMN']);
write2file([Path,'\src\fem\mex\localsNH']);
%write2file(mexPath);
else
addPath([Path,'\scripts\fem\']);
addPath([Path,'\scripts\fem\2D']);
addPath([Path,'\scripts\fem\3D']);
addPath([Path,'\scripts\opt']);
addPath([Path,'\src\fem']);
addPath([Path,'\src\fem\tools\']);
addPath([Path,'\src\fem\tools\tpswarp']);
addPath([Path,'\src\fem\tools\mma']);
addPath([Path,'\src\fem\mex\']);
addPath([Path,'\src\fem\materials\']);
addPath([Path,'\src\fem\materials\samples']);
addPath([Path,'\src\fem\mex\nonlinearstrainmat']);
addPath([Path,'\src\fem\mex\polardecomposition']);
addPath([Path,'\src\fem\mex\localsYH']);
addPath([Path,'\src\fem\mex\localsMN']);
addPath([Path,'\src\fem\mex\localsNH']);
pause(.3);

x = checkLibary('fem.lib',@(x) femPathConfirm);

if x
    fprintf(FIT,'%% fem.lib installed succesfully! \n');
else
    fprintf(FIT,'%% fem.lib missing! \n');
end

end
end
% ---------------------------------------------------------------- DYNAMICS
function x = includeDynamicModel(Path,Request)
global FID FIT

if Request == 1
fprintf(FID,'%% mdl.lib \n');
write2file([Path,'\scripts\mdl']);
write2file([Path,'\livescripts\mdl']);
write2file([Path,'\src\model']);
write2file([Path,'\src\model\tools']);
write2file([Path,'\src\model\tools\liegroup']);
write2file([Path,'\src\model\tools\shapefnc']);
write2file([Path,'\src\model\mex\forwardkin']);
write2file([Path,'\src\model\mex\lagrangian']);
%write2file([Path,'\src\model\tools\solver']);
% mexPath = genpath([Path,'\src\model\mex']);
% write2file(mexPath);
else
addPath([Path,'\scripts\mdl']);
addPath([Path,'\livescripts\mdl']);
addPath([Path,'\src\model']);
addPath([Path,'\src\model\tools']);
addPath([Path,'\src\model\tools\liegroup']);
addPath([Path,'\src\model\tools\shapefnc']);
addPath([Path,'\src\model\tools']);
addPath([Path,'\src\model\mex\forwardkin']);
addPath([Path,'\src\model\mex\lagrangian']);
% mexPath = genpath([Path,'/src/model/mex']);
% addPath(mexPath);
pause(.3);

x = checkLibary('model.lib',@(x) modelPathConfirm);

if x, fprintf(FIT,'%% model.lib installed succesfully! \n');
else, fprintf(FIT,'%% model.lib missing! \n');
end

end
end
% ----------------------------------------------------------------- CONTROL
function x = includeControl(Path,Request)
global FID FIT

if Request == 1
fprintf(FID,'%% bdog.lib \n');
write2file([Path,'\scripts\bdog']);
write2file([Path,'\src\bdog']);
write2file([Path,'\src\bdog\tools']);
write2file([Path,'\src\vision']);
write2file([Path,'\src\vision\tools']);
else
addPath([Path,'\scripts\bdog']);
addPath([Path,'\src\bdog']);
addPath([Path,'\src\bdog\tools']);
addPath([Path,'\src\vision']);
addPath([Path,'\src\vision\tools']);
pause(.3);

x = checkLibary('control.lib',@(x) bdogPathConfirm);

if x
    fprintf(FIT,'%% control.lib installed succesfully! \n');
else
    fprintf(FIT,'%% control.lib missing! \n');
end

end

end

