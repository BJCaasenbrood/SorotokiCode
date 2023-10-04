function [name, date] = getPackagesInformationMPM()
    mpmPath = which('mpm');
    mpmPackageFolder = [mpmPath(1:end-6),'/mpm-packages/'];
    try
        MPM = load([mpmPackageFolder,'mpm.mat']);
        name = {MPM.packages.name};
        date = {MPM.packages.downloadDate};
    catch
        name = [];
        date = [];
    end
end
