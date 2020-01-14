function repo(Request)

if nargin < 1, Request = ''; end

switch(Request)
    case('demo');      DirectoryDemos;
    case('demos');     DirectoryDemos;
    case('model');     DirectoryModels;
    case('models');    DirectoryModels;
    case('stl');       DirectoryModels;
    case('mat');       DirectoryMaterials;
    case('material');  DirectoryMaterials;
    case('materials'); DirectoryMaterials;
    case('sdf');       DirectorySDF;
    case('imp');       DirectorySDF;
    case('implicit');  DirectorySDF;
    case('color');     DirectoryColors;
    case('colors');    DirectoryColors;
    case('colormap');  DirectoryColormap;
    case('colormaps'); DirectoryColormap;
    otherwise;         BaseRepo;
end
end

function BaseRepo
cout('Data:')
dir('./data/')
cout('Demos:')
dir('./scripts/')
end

function DirectoryColors
dir('.\data\colors\**\*.m');
end

function DirectorySDF
dir('.\examples\resources\implicit\**\*.m');
end

function DirectoryModels
dir('.\data\stl\**\*.stl');
end

function DirectoryMaterials
dir('.\data\matcap\**\*.m');
end


function DirectoryDemos
dir('.\examples\demos\**\*.m');
end