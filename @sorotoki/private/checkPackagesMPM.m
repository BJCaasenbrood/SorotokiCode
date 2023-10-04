function checkPackagesMPM(reqPackages,mpmPackages)

    global log auto_approve;
    log.info('Checking for updates: MPM package library'); 
    tf = ismember(reqPackages,mpmPackages);

    if all(tf)
        log.info(msg('mpm_complete'));
    else
        log.info(msg('mpm_incomplete'));
        I = find(~tf).';
        log.list('',reqPackages(I(:)));

        if ~auto_approve
            log.help('Install missing MPM packages?');
            reply = input(i18n('confirm'), 's');
            if isempty(reply)
                reply = i18n('confirm_yes');
            end
            if ~strcmpi(reply(1), i18n('confirm_yes'))
                log.info(i18n('confirm_nvm'));
                return;
            end
        end
    
        misPackages = reqPackages(find(~tf).');
        log.info('Calling MPM installer -- installing req. packages');
        log.hline();
        log.setHide(true);
        for i = 1:numel(misPackages)
            installMissingPackageMPM(misPackages{i}); 
            log.setHide(false);
            log.hline();
            log.setHide(true);
        end
    end
end