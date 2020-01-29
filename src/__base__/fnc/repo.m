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

fprintf('=====================================================\n')
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
    otherwise;         BaseRepo;
end
fprintf('=====================================================\n')
end

function BaseRepo
dir('.\data\color\')
dir('.\data\colormap\')
dir('.\data\matcap\')
dir('.\scripts\**')
end

function DirectoryColors
dir('.\data\colors\**\*.m');
end

function DirectorySDF
dir('.\examples\resources\implicit\**\*.m');
end

function DirectoryModels
dir('.\data\stl\**\*.stl');
dir('.\data\stl\**\*.obj');
end

function DirectoryMaterials
dir('.\data\matcap\**\*.m');
end

function DirectoryDemos
dir('.\scripts\**\*.m');
end