function installSorotoki(reqPackages,soroPackages)
    global log;

    installDir = pwd;
    mpmPath    = which('mpm');
    [mpmPackages, mpmPackagedate] = getPackagesInformationMPM();

    if ~isfile('sorotoki.log')
        log.debug('Deleting sorotoki.log');
        delete sorotoki.log

        log.debug('Opening diary sorotoki.log');
        diary sorotoki.log
        log.info(['Installation dir: ', installDir]);
        log.info(['Directory mpm found: ', mpmPath(1:end-6)]);
    end
    
    log.debug('Checking for packages packages');
    checkPackagesMPM(reqPackages,mpmPackages);

    log.debug('Checking for sorotoki packages');
    checkSoroPackages(soroPackages,mpmPackages);
    
    diary off;
end