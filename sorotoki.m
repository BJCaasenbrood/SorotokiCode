function sorotoki(arg)
% -------------------------------------------------------------------------
if nargin < 1, clc; clear; arg = 'install'; end
% ------------------------------------------------------------------------

switch(arg)
    case('install');  setupToolkit;
    case('demo');     showDemo;
    case('cd');       getPath;
end

end

% -------------------------------------------------------------- ESSENTIALS
function showDemo
clr;
set = {'Meshing a 2D-circle',...
       'Rendering stanford bunny',...
       'Nonlinear finite element',...
       '3D nonlinear finite element',...
       'Topology optimization of Pneu-net'};

request = action(set);

switch(request)
    case(1); open mesh_circle;
    case(2); open preview_bunny;
    case(3); open fem_bellow;
    case(4); open fem_twisting_beam;
    case(5); open opt_pneunet;
    otherwise; warning('Please select a demo from the list above.');
end

end

% -------------------------------------------------------------- ESSENTIALS
function setupToolkit  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all; warning off; beep off;

addpath('src');
addpath('src/__version__');
addpath('src/__base__');
addpath('src/__base__/fnc');

skipUpdate = false;
Path = cd;
DisplayLogo;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(['* Thank you for installing SOROTOKI, we will',...
    ' start the installation shortly \n']); pause(1.0);

% if exist('BuildVersion.m','file')
%      % File exists.
% else
%      % File does not exist.
% end

fprintf('* Starting installation - SOROTOKI');
verFolder = 'config';
if ~exist(verFolder, 'dir')
mkdir(verFolder);
fprintf(['* Created directory ',verFolder, ' \n']);
fprintf(['* Directory ',verFolder, ' added to path \n']);
else
fprintf(['* Directory ',verFolder, ' added to path \n']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AddPath([Path,verFolder]);

if exist([Path,'/config/vernum.m'], 'file')
    delete([Path,'/config/vernum.m']); 
end

if ~pingserver
    cout('err','No internet connection! '); 
    str = action({'(y)es, continue without connection',...
            '(n)o, stop installation'},'s');
    if strcmp(str,'N'), return; end
    skipUpdate = true;
end

if ~skipUpdate
fprintf('* Getting version_file_check from Git repository config/vernum.m \n');
url = ['https://raw.githubusercontent.com/BJCaasenbrood/',...
    'SorotokiCode/master/config/vernum.m'];
filename = [verFolder,'/vernum.m'];
websave(filename,url);

fprintf(['* Succesfully downloaded contents', filename, '\n']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
libs(1) = IncludeBase(Path,0);
libs(2) = IncludeGraphicsModel(Path,0);
libs(3) = IncludeMesh(Path,0);
libs(4) = IncludeFiniteElement(Path,0);
libs(5) = IncludeDynamicModel(Path,0);
libs(6) = IncludeControl(Path,0);

if min(libs) == 1
fprintf('\n* Libary check completed - all libaries are up-to-date - \n');
else
cout('err','\n* Libary check completed - some libaries are outdated!\n');
fprintf('* Get the latest version at: ');

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
    global FID;
    StartUpFile = [userpath,'/startup.m'];
    delete(StartUpFile);
    FID = fopen(StartUpFile,'w');
    
    fprintf('* Assigning fixed search paths to MATLAB: \n');
    if libs(1), IncludeBase(Path,1); end
    if libs(2), IncludeGraphicsModel(Path,1); end
    if libs(3), IncludeMesh(Path,1); end
    if libs(4), IncludeFiniteElement(Path,1); end
    if libs(5), IncludeDynamicModel(Path,1); end
    if libs(6), IncludeControl(Path,1); end
    
    Request = input(['\n* Do you want to setup the SorotokiCode folder',...
        ' as the main directory (y/n)'],'s');

    switch(Request)
        case('y'), fprintf(FID,'cdsoro; \n');
    end
    
    fclose('all');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = ver;
ToolboxesInstalled = {S.Name}';
Requisites = {'Image Processing Toolbox';'Partial Differential Equation Toolbox'};
HTMLLink = {'<a href="https://nl.mathworks.com/products/image.html">click here to install</a> \n',...
            '<a href="https://nl.mathworks.com/products/pde.html">click here to install</a> \n'};
memb = ismember(Requisites,ToolboxesInstalled);
cout('Text','* Checking for prerequisit toolboxes...\n'); pause(1);
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

cout('Text','* Checking for Cmake...\n'); pause(1);

[status,~] = system("cmake --version");
if  status~= 0
    cout('Error', '* Missing Cmake! \n');    
else
    cout('green','* Cmake compiler founded\n')
end
pause(0.75);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cout('Text', '\n');
cout('green','* INSTALLATION DONE! \n');
cout('Text', '\n');
cout('* SOROTOKI toolkit is succesfully installed and ready-to-use! \n');

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

% -------------------------------------------------------------- ESSENTIALS
function Path = getPath(print)
if nargin < 1; print = false;end
Path = cd;
if print, disp(['* cd: ',Path]); end
end

% -------------------------------------------------------------- ESSENTIALS
function x = IncludeBase(Path,Request)
global FID
cout(['* Adding SOROTOKI libraries to path, ',...
    'this might take a minute...\n\n']);

if Request == 1
%str = strrep(Path,'/','\');
fprintf(FID,['%%!INSTALDIR:',strrep(Path,'\','/'),' \n']);
fprintf(FID,'%% base.lib \n');
%WriteToFile(Path);
WriteToFile([Path,'\config\']);
%WriteToFile([Path,'\src\']);
WriteToFile([Path,'\src\__version__']);
WriteToFile([Path,'\src\__base__']);
WriteToFile([Path,'\src\__base__\fnc']);
WriteToFile([Path,'\data\']);
WriteToFile([Path,'\data\color']);
WriteToFile([Path,'\data\colormap']);
WriteToFile([Path,'\data\matcap']);
WriteToFile([Path,'\data\matcap\img']);
WriteToFile([Path,'\data\stl']);
WriteToFile(genpath([Path,'/scripts']));
else
AddPath(Path);
AddPath([Path,'\config\']);
AddPath([Path,'\src\__version__']);
AddPath([Path,'\src\__base__']);
AddPath([Path,'\src\__base__\fnc']);
AddPath([Path,'\data\color']);
AddPath([Path,'\data\colormap']);
AddPath([Path,'\data\matcap']);
AddPath([Path,'\data\matcap\img']);
AddPath([Path,'\data\stl']);
AddPath(genpath([Path,'/scripts']));
pause(.3);
x = basePathConfirm;
end
end

% ---------------------------------------------------------------- GRAPHICS
function x = IncludeGraphicsModel(Path,Request)
global FID

if Request == 1
fprintf(FID,'%% gmodel.lib \n');
WriteToFile([Path,'\src\gmodel']);
WriteToFile([Path,'\src\gmodel\tools\']);
WriteToFile([Path,'\src\gmodel\matcap\']);
WriteToFile([Path,'\src\gmodel\matcap\img']);
else
AddPath([Path,'\src\gmodel']);
AddPath([Path,'\src\gmodel\tools\']);
AddPath([Path,'\src\gmodel\matcap\']);
AddPath([Path,'\src\gmodel\matcap\img\']);
AddPath([Path,'\src\gmodel\matcap\tools\']);
pause(.3);
x = graphicsmodelPathConfirm;
end
end

% -------------------------------------------------------------------- MESH
function x = IncludeMesh(Path,Request)
global FID

if Request == 1
fprintf(FID,'%% mesh.lib \n');
WriteToFile([Path,'\src\mesh']);
WriteToFile([Path,'\src\mesh\tools\']);
WriteToFile([Path,'\src\mesh\shapes\']);
WriteToFile([Path,'\src\mesh\operators\']);
else
AddPath([Path,'\src\mesh']);
AddPath([Path,'\src\mesh\tools\']);
AddPath([Path,'\src\mesh\shapes\']);
AddPath([Path,'\src\mesh\operators\']);
pause(.3);
x = meshPathConfirm;
end
end

% --------------------------------------------------------------------- FEM
function x = IncludeFiniteElement(Path,Request)
global FID

if Request == 1
fprintf(FID,'%% fem.lib \n');
WriteToFile([Path,'\src\fem']);
WriteToFile([Path,'\src\fem\tools\']);
WriteToFile([Path,'\src\fem\tools\tpswarp']);
WriteToFile([Path,'\src\fem\tools\mma']);
WriteToFile([Path,'\src\fem\materials\']);
WriteToFile([Path,'\src\fem\materials\samples']);
else
AddPath([Path,'\src\fem']);
AddPath([Path,'\src\fem\tools\']);
AddPath([Path,'\src\fem\tools\tpswarp']);
AddPath([Path,'\src\fem\tools\mma']);
AddPath([Path,'\src\fem\materials\']);
AddPath([Path,'\src\fem\materials\samples']);
pause(.3);
x = femPathConfirm;
end
end

% ---------------------------------------------------------------- DYNAMICS
function x = IncludeDynamicModel(Path,Request)
global FID

if Request == 1
fprintf(FID,'%% mdl.lib \n');
WriteToFile([Path,'\src\model']);
WriteToFile([Path,'\src\model\tools']);
WriteToFile([Path,'\src\model\tools\solver']);
else
AddPath([Path,'\src\model']);
AddPath([Path,'\src\model\tools']);
pause(.3);
x = modelPathConfirm;
end
end

% ----------------------------------------------------------------- CONTROL
function x = IncludeControl(Path,Request)
global FID

if Request == 1
fprintf(FID,'%% bdog.lib \n');
WriteToFile([Path,'\src\bdog']);
WriteToFile([Path,'\src\bdog\tools']);
else
AddPath([Path,'\src\bdog']);
AddPath([Path,'\src\bdog\tools']);
pause(.3);
x = bdogPathConfirm;
end
end

