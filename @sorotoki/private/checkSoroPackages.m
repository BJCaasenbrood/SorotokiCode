function checkSoroPackages(soroPackages,mpmPackages)

    global log

    log.info('Checking for updates: Sorotoki library'); 
    tf = ismember(soroPackages,mpmPackages);
    
    if all(tf)
        log.info(msg('soro_complete'));
    else
        log.info(msg('soro_incomplete'));
        I = find(~tf).';
        log.list('',soroPackages(I(:)));

        log.help('Install missing MPM packages?');
        reply = input(i18n('confirm'), 's');
        if isempty(reply)
            reply = i18n('confirm_yes');
        end
        if ~strcmpi(reply(1), i18n('confirm_yes'))
            log.info(i18n('confirm_nvm'));
            return;
        end
    
        misPackages = soroPackages(find(~tf).');
        log.info('Calling MPM installer -- installing Sorotoki packages');
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
