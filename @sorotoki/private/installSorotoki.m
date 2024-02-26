function installSorotoki(reqPackages,soroPackages)
    % global log;

    installDir = pwd;

    mpiPath = which('mpi');
    [mpiPackages, ~] = getPackagesInformationMPI();

    if ~isfile('sorotoki.log')
        % log.debug('Deleting sorotoki.log');
        delete sorotoki.log

        % log.debug('Opening diary sorotoki.log');
        diary sorotoki.log
        disp(['Installation dir: ', installDir]);
        disp(['Directory mpm found: ', mpiPath(1:end-6)]);
    end
    
    % log.debug('Checking for packages packages');
    checkPackagesMPI(reqPackages,mpiPackages);

    % log.debug('Checking for sorotoki packages');
    checkSoroPackages(soroPackages,mpiPackages);
    diary off;
end