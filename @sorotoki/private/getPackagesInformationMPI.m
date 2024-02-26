function [name, date] = getPackagesInformationMPI()
    mpiPath = which('mpi');
    mpiPackageFolder = [mpiPath(1:end-6),'/mpi-packages/'];
    try
        MPI = load([mpiPackageFolder,'mpi.mat']);
        name = {MPI.packages.name};
        date = {MPI.packages.downloadDate};
    catch
        name = {''};
        date = {''};
    end
end
