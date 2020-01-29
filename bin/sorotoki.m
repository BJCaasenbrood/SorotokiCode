function sorotoki(arg)
% -------------------------------------------------------------------------
if nargin < 1, clc; clear; arg = 'demo'; end
if ~exist('config/__init__.txt','file'), arg = 'install'; end
% ------------------------------------------------------------------------

switch(arg)
    case('install');  setupToolkit;
    case('demo');     showDemo;
end

end

function showDemo
ccc;
set = {'bunny preview demo',...
       'dragon preview demo',...
       'advanced FV-mesh preview demo',...
       'implicit meshing demo',...
       'linear fem demo',...
       'nonlinear fem demo',...
       'topology design demo'};

request = CallAction(set);
switch(request)
    case(1); preview_bunny;
    case(2); preview_dragon;
    case(3); preview_demo;
    case(4); mesh_demo;
    case(5); fem_demo;
    case(6); nlfem_demo;
    case(7); design_demo;
    otherwise; CallWarning('please select demo from the list above.');
end

end

function setupToolkit  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all; warning off; beep off;

addpath('src');
addpath('src/__version__');
addpath('src/__base__');

Path = getPath;
DisplayLogo;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(['* Thank you for installing SOROTOKI, we will',...
    ' start the installation shortly \n']); pause(1.0);

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
addpath([Path,verFolder]);

if exist([Path,'/config/vernum.m'], 'file')
    delete([Path,'/config/vernum.m']); 
end

if ~PingInternet, error('Enable your internet connection!'); end

fprintf('* getting version_file_check from Git repository config/vernum.m \n');
url = ['https://raw.githubusercontent.com/BJCaasenbrood/',...
    'SorotokiCode/master/config/vernum.m'];
filename = [verFolder,'/vernum.m'];
websave(filename,url);

fprintf(['* Succesfully downloaded contents', filename, '\n']);
%fprintf(['* Creating installation log file', filename, '\n']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
libs(1) = IncludeBase(Path,0);
libs(2) = IncludeGraphicsModel(Path,0);
libs(3) = IncludeMesh(Path,0);
libs(4) = IncludeFiniteElement(Path,0);
libs(5) = IncludeDynamicModel(Path,0);

if min(libs) == 1
fprintf('\n* Libary check completed - all libaries are up-to-date \n');
else
cprintf('err','\n* Libary check completed - some libaries are outdated!\n');
fprintf('* Get the latest version at: ');

fprintf(['<a href="https://github.com/BJCaasenbrood/SorotokiCode">',...
    'https://github.com/BJCaasenbrood/Sorotoki</a> \n'])

Request = input('  Do you want to continue the installation? (y/n) ','s');

switch(Request)
    case('y'); fprintf('* Proceeding toolkit installation... \n'); pause(1);
    case('n'); error('installation aborted!');
    otherwise; error('installation aborted!');
end

end

Request = input(' Do you want to create a startup file? (y/n) ','s');
bool = 1;
switch(Request)
    case('y'); fprintf('* Proceeding generation startup.m file... \n'); pause(1);
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
    
    fclose('all');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
startup;
cout('Text', '\n');
cout('green','* INSTALLATION DONE! \n');
cout('Text', '\n');
cout('SOROTOKI toolkit is succesfully installed and ready-to-use! ');

pause(.01);
cout(['The documentation can be found in doc/SorotokiManual.pdf',...
    '. For more in-\nformation on the Soft Robotics Toolkit, visit',...
    ' the GitHub repository at: \n\n']) 
cout('Text',['\t <a href="https://github.com/BJCaasenbrood/Sorotoki">',...
    'https://github.com/BJCaasenbrood/SorotokiCode</a>\n']);
cout('Text', '\n');
cout('Text', '* To get started with demonstrations, type ');
cout('String', 'sorotoki ');
cout('Text', 'or '); 
cout('String', 'sorotoki(''demo'') \n');
warning on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function Path = getPath(print)
if nargin < 1; print = false;end
Path = cd;
if print, disp(['* cd: ',Path]); end
end

% -------------------------------------------------------------- ESSENTIALS
function x = IncludeBase(Path,Request)
global FID

cout(['* Adding SOROTOKI libraries to path, ',...
    'this might take a minute...\n']);

if Request == 1
fprintf(FID,'%% base.lib \n');
WriteToFile([Path,'\config\']);
WriteToFile([Path,'\src\__version__']);
WriteToFile([Path,'\src\__base__']);
WriteToFile([Path,'\src\__base__\fnc']);
WriteToFile([Path,'\scripts\']);
WriteToFile([Path,'\data\color']);
WriteToFile([Path,'\data\colormap']);
WriteToFile([Path,'\data\matcap']);
WriteToFile([Path,'\data\matcap\img']);
WriteToFile([Path,'\data\stl']);
else
addpath([Path,'\config\']);
addpath([Path,'\src\__version__']);
addpath([Path,'\src\__base__']);
addpath([Path,'\src\__base__\fnc']);
addpath([Path,'\scripts\']);
addpath([Path,'\data\color']);
addpath([Path,'\data\colormap']);
addpath([Path,'\data\matcap']);
addpath([Path,'\data\matcap\img']);
addpath([Path,'\data\stl']);
pause(.3);
x = basePathConfirm;
end
end

% -------------------------------------------------------------- ESSENTIALS
function x = IncludeGraphicsModel(Path,Request)
global FID

if Request == 1
fprintf(FID,'%% gmodel.lib \n');
WriteToFile([Path,'\src\gmodel']);
WriteToFile([Path,'\src\gmodel\tools\']);
else
addpath([Path,'\src\gmodel']);
addpath([Path,'\src\gmodel\tools\']);
pause(.3);
x = graphicsmodelPathConfirm;
end
end

% -------------------------------------------------------------- ESSENTIALS
function x = IncludeMesh(Path,Request)
global FID

if Request == 1
fprintf(FID,'%% mesh.lib \n');
WriteToFile([Path,'\src\mesh']);
WriteToFile([Path,'\src\mesh\tools\']);
WriteToFile([Path,'\src\mesh\shapes\']);
WriteToFile([Path,'\src\mesh\operators\']);
else
addpath([Path,'\src\mesh']);
addpath([Path,'\src\mesh\tools\']);
addpath([Path,'\src\mesh\shapes\']);
addpath([Path,'\src\mesh\operators\']);
pause(.3);
x = meshPathConfirm;
end
end

% -------------------------------------------------------------- ESSENTIALS
function x = IncludeFiniteElement(Path,Request)
global FID

if Request == 1
fprintf(FID,'%% fem.lib \n');
WriteToFile([Path,'\src\fem']);
WriteToFile([Path,'\src\fem\tools\']);
WriteToFile([Path,'\src\fem\tools\mma']);
WriteToFile([Path,'\src\fem\materials\']);
else
addpath([Path,'\src\fem']);
addpath([Path,'\src\fem\tools\']);
addpath([Path,'\src\fem\tools\mma']);
addpath([Path,'\src\fem\materials\']);
pause(.3);
x = femPathConfirm;
end
end

% -------------------------------------------------------------- ESSENTIALS
function x = IncludeDynamicModel(Path,Request)
global FID

if Request == 1
fprintf(FID,'%% mdl.lib \n');
WriteToFile([Path,'\src\model']);
WriteToFile([Path,'\src\model\tools']);
else
addpath([Path,'\src\model']);
addpath([Path,'\src\model\tools']);
pause(.3);
x = modelPathConfirm;
end
end

