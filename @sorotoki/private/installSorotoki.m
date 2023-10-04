function installSorotoki(reqPackages,soroPackages)
    % global log;

    installDir = pwd;
    mpmPath    = which('mpm');
    [mpmPackages, ~] = getPackagesInformationMPM();

    if ~isfile('sorotoki.log')
        % log.debug('Deleting sorotoki.log');
        delete sorotoki.log

        % log.debug('Opening diary sorotoki.log');
        diary sorotoki.log
        disp(['Installation dir: ', installDir]);
        disp(['Directory mpm found: ', mpmPath(1:end-6)]);
    end
    
    % log.debug('Checking for packages packages');
    checkPackagesMPM(reqPackages,mpmPackages);

    % log.debug('Checking for sorotoki packages');
    checkSoroPackages(soroPackages,mpmPackages);
    diary off;
end