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
DisplayLogo

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

fprintf('* getting version_file_check from repository config/vernum.m \n');
url = ['https://raw.githubusercontent.com/BJCaasenbrood/',...
    'SorotokiCode/master/config/vernum.m'];
filename = [verFolder,'/vernum.m'];
websave(filename,url);

fprintf(['* Succesfully downloaded contents', filename, '\n']);
%fprintf(['* Creating installation log file', filename, '\n']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
libs(1) = IncludeBase(Path,0);
libs(2) = IncludeObject(Path,0);

if min(libs) == 1
fprintf('\n* Libary check completed - all libaries are up-to-date \n');
else
cprintf('err','\n* Libary check completed - some libaries are outdated!\n');
fprintf('* Get the latest version at: ');

fprintf(['<a href="https://github.com/BJCaasenbrood/SorotokiCode">',...
    'https://github.com/BJCaasenbrood/Sorotoki</a> \n'])

Request = input('  Do you want to continue the installation? (y/n) ','s');

switch(Request)
    case('y'); fprintf('* Proceeding toolkit installation... \n');pause(1);
    case('n'); error('installation aborted!');
    otherwise; error('installation aborted!');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global FID;
StartUpFile = [userpath,'/startup.m'];
delete(StartUpFile);
FID = fopen(StartUpFile,'w');

fprintf('* Assigning fixed search paths to MATLAB: \n');
if libs(1), IncludeEssentials(Path,1); end
if libs(2), IncludeBase(Path,1); end

fclose('all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
startup;

cout('green','* INSTALLATION DONE! \n');
cout('Text', '\n');
cout('SOROTOKI toolkit is succesfully installed and ready-to-use!');

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

cout(['* Adding SOROTOKI libraries to path,',...
    'this might take a minute...\n']);

if Request == 1
fprintf(FID,'%% base.lib \n');
WriteToFile([Path,'/config/']);
WriteToFile([Path,'/src/__version__']);
WriteToFile([Path,'/src/__base__']);
WriteToFile([Path,'/data/colors']);
WriteToFile([Path,'/data/matcap']);
WriteToFile([Path,'/data/matcap/img']);
else
addpath([Path,'/config/']);
addpath([Path,'/src/__version__']);
addpath([Path,'/src/__base__']);
addpath([Path,'/data/colors']);
addpath([Path,'/data/matcap']);
addpath([Path,'/data/matcap/img']);
pause(.3);
x = basePathConfirm;
end
end

% -------------------------------------------------------------- ESSENTIALS
function x = IncludeObject(Path,Request)
global FID

if Request == 1
fprintf(FID,'%% object.lib \n');
WriteToFile([Path,'/src/object']);
WriteToFile([Path,'/src/object/tools']);
else
addpath([Path,'/src/object']);
addpath([Path,'/src/object/tools']);
pause(.3);
x = objectPathConfirm;
end
end

% -------------------------------------------------------------- ESSENTIALS
function x = IncludeEssentials(Path,Request)
global FID
if Request == 1
fprintf(FID,'%% output.lib \n');
WriteToFile([Path,'\src\output']);
else
addpath([Path,'\src\output']);
pause(.3);
x = interfacePathConfirm;
end
end

% ----------------------------------------------------------------- MESHING
function x = IncludePreview(Path,Request)
global FID;
if Request == 1
fprintf(FID,'%% preview.lib \n');
WriteToFile([Path,'\sorotoki_preview_master']);
WriteToFile([Path,'\sorotoki_preview_master\tools']);
WriteToFile([Path,'\sorotoki_preview_master\tools\colors']);
WriteToFile([Path,'\examples\resources\models'])
WriteToFile([Path,'\examples\resources\cubemaps'])
WriteToFile([Path,'\examples\resources\matcaps'])
x = 1;
else
addpath([Path,'\sorotoki_preview_master']);
x = previewPathConfirm;
end
end

% ---------------------------------------------------------- PHYSICS ENGINE
function x = IncludeBlender(Path,Request)
global FID;
if Request == 1
fprintf(FID,'%% blender.lib \n');
WriteToFile([Path,'\sorotoki_blender_master']);
x = 1;
else
addpath([Path,'\sorotoki_blender_master']);
x = blenderPathConfirm;
end
end

% ----------------------------------------------------------------- MESHING
function x = IncludeMesher(Path,Request)
global FID;
if Request == 1
fprintf(FID,'%% mesher.lib \n');
WriteToFile([Path,'\sorotoki_mesher_master']);
WriteToFile([Path,'\sorotoki_mesher_master\tools']);
WriteToFile([Path,'\sorotoki_mesher_master\tools\shapes']);
WriteToFile([Path,'\sorotoki_mesher_master\tools\operators']);
WriteToFile([Path,'\examples\resources\implicit'])
else
addpath([Path,'\sorotoki_mesher_master']);
x = mesherPathConfirm;
end
end

% --------------------------------------------------------- FINITE ELEMENTS
function x = IncludeFiniteElements(Path,Request)
global FID;
if Request == 1
fprintf(FID,'%% fem.lib \n');
WriteToFile([Path,'\sorotoki_finite_elements_master']);
WriteToFile([Path,'\sorotoki_finite_elements_master\tools']);
WriteToFile([Path,'\sorotoki_finite_elements_master\materials']);
else
addpath([Path,'\sorotoki_finite_elements_master']);
x = fePathConfirm;
end
end

% --------------------------------------------------- TOPOLOGY OPTIMIZATION
function x = IncludeTopologyOptimzation(Path,Request)
global FID;
if Request == 1
fprintf(FID,'%% topology.lib \n');   
WriteToFile([Path,'\sorotoki_topology_master']);
WriteToFile([Path,'\sorotoki_topology_master\tools']);
else
addpath([Path,'\sorotoki_topology_master']);
x = topologyOptPathConfirm;
end
end

% ------------------------------------------------------------------ FORMER
function x = IncludeFormer(Path,Request)
global FID;
if Request == 1
fprintf(FID,'%% former.lib \n');
WriteToFile([Path,'\sorotoki_former_master']);
WriteToFile([Path,'\sorotoki_former_master\tools']);
else
addpath([Path,'\sorotoki_former_master']);
x = formerPathConfirm;
end
end

% -------------------------------------------------------- BALLOONDOG BOARD
function x = IncludeBalloonDog(Path,Request)
global FID;
if Request == 1
fprintf(FID,'%% balloondog.lib \n');
WriteToFile([Path,'\sorotoki_balloondogPi_master']);
WriteToFile([Path,'\sorotoki_balloondogPi_master\SSH']);
else
addpath([Path,'\sorotoki_balloondogPi_master']);
addpath([Path,'\sorotoki_balloondogPi_master\SSH']);
x = balloondogPathConfirm;
end
end

% -------------------------------------------------- ADDITIVE MANUFACTURING
function x = IncludePrinter(Path)
addpath([Path,'\sorotoki_printer_master']);
printerPathConfirm;
end

% ----------------------------------------------------------- MAGNETICS FEM
function x= IncludeMagnetics(Path,Request)
global FID;
if Request == 1
fprintf(FID,'%% balloondog.lib \n');
WriteToFile([Path,'\sorotoki_magnetics_master']);
WriteToFile([Path,'\sorotoki_magnetics_master\tools']);
WriteToFile([Path,'\src\examples\domains\magnet']);
else
addpath([Path,'\sorotoki_magnetics_master']);
addpath([Path,'\sorotoki_magnetics_master\tools']);
addpath([Path,'\src\examples\domains\magnet']);
x = magneticsPathConfirm;
end
end

% ----------------------------------------------------------- MAGNETICS FEM
function x= IncludeClasses(Path,Request)
global FID;
if Request == 1
fprintf(FID,'%% balloondog.lib \n');
WriteToFile([Path,'\sorotoki_classes_master']);
else
addpath([Path,'\sorotoki_classes_master']);
x = classesPathConfirm;
end
end


