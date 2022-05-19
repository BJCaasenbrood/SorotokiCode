function repo(Request)
%REPO Repository of the soft robotics toolkit
%   REPO(X) gives an overview of all files that are
%   contained in the (X) map of SOROTOKI toolkit.
%
%   REPO('demo') returns all the demo files,
%   other commands are given as:
%
%   Class support for inputs:
%      'demo', 'demos', 'scripts'
%      'model', 'models', 'stl', 'obj',
%      'mat', 'matcap', 'material', 
%      'sdf', 'implicit', 'function',
%      'color', 'colors', 
%      'colormap', 'colormaps'.
%
%   See also CD, CDSORO, PING
%
%   For more information: 
%       https://github.com/BJCaasenbrood/SorotokiCode

%   Copyright 2018-2023 B.J.Caasenbrood 

if nargin < 1, Request = ''; end

switch(Request)
    case('demo');      DirectoryDemos;
    case('demos');     DirectoryDemos;
    case('model');     DirectoryModels;
    case('models');    DirectoryModels;
    case('stl');       DirectoryModels;
    case('obj');       DirectoryModels;
    case('mat');       DirectoryMaterials;
    case('material');  DirectoryMaterials;
    case('matcap');    DirectoryMaterials;
    case('sdf');       DirectorySDF;
    case('function');  DirectorySDF;
    case('functions'); DirectorySDF;
    case('implicit');  DirectorySDF;
    case('color');     DirectoryColors;
    case('colors');    DirectoryColors;
    case('colormap');  DirectoryColormap;
    case('colormaps'); DirectoryColormap;
    case('opt');       DirectoryOptimization;
    otherwise;         BaseRepo;
end

end

function BaseRepo
p = cdsoro;
dir([p,'\data\color\']);
dir([p,'\data\colormap\']);
dir([p,'\data\matcap\']);
dir([p,'\scripts\**']);
end

function DirectoryColors
p = cdsoro;
dir([p,'\data\colors\**\*.m']);
end

function DirectorySDF
p = cdsoro;
dir([p,'\examples\resources\implicit\**\*.m']);
end

function DirectoryModels
p = cdsoro;
dir([p,'\data\stl\**\*.stl']);
dir([p,'\data\stl\**\*.obj']);
end

function DirectoryMaterials
p = cdsoro;
dir([p,'/data/matcap/**/*.m']);
end

function DirectoryDemos
p = cdsoro;
dir([p,'\scripts\**\*.m']);
end

function DirectoryOptimization
p = cdsoro;
dir([p,'\scripts\opt\*.m']);
end

function DirectoryColormap
p = cdsoro;
dir([p,'\data\colormap\*.m']);
end