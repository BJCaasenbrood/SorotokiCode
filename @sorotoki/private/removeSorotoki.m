function removeSorotoki(soroPackages)
    
    global log;
    
    warning off;

    log.info('Calling MPM installer -- removing Sorotoki');
    log.hline();
    for i = 1:numel(soroPackages)
        installMissingPackageMPM(soroPackages{i},'uninstall'); 
    end
    log.hline();

    try
        rmdir('lib');
        rmdir('assets');
        delete sorotoki.log
    end

    warning on;
end