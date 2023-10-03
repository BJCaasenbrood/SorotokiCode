function gotoSorotokiFolder
    log.debug('Finding installation log file: sorotoki.log');
    installPath = whereis('sorotoki.log');
    log.debug(['Found: ',installPath,'/sorotoki.log']);
    log.info(['cd into ', installPath(1:end-10)]);
    cd(installPath(1:end-10));
end