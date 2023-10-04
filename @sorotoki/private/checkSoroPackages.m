function checkSoroPackages(soroPackages,mpmPackages)

    global auto_approve

    disp('Checking for updates: Sorotoki library'); 
    tf = ismember(soroPackages,mpmPackages);
    
    if all(tf)
        disp(msg('soro_complete'));
    else
        disp(msg('soro_incomplete'));
        I = find(~tf).';
        for i = 1:numel(I)
            disp(['- ', soroPackages{I(i)}]);
        end

        if ~auto_approve
            disp('Install missing MPM packages?');
            reply = input(i18n('confirm'), 's');
            if isempty(reply)
                reply = i18n('confirm_yes');
            end
            if ~strcmpi(reply(1), i18n('confirm_yes'))
                log.info(i18n('confirm_nvm'));
                return;
            end
        end
    
        lineStr = repmat('━', 1, 40);
        misPackages = soroPackages(find(~tf).');
        disp('Calling MPM installer -- installing Sorotoki packages');
        disp(lineStr);
        for i = 1:numel(misPackages)
            installMissingPackageMPM(misPackages{i}); 
            disp(lineStr);
        end
    end
end
