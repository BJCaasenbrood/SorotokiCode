function checkPackagesMPM(reqPackages,mpmPackages)

    global auto_approve;
    disp('Checking for updates: MPM package library'); 
    tf = ismember(reqPackages,mpmPackages);

    if all(tf)
        disp(msg('mpm_complete'));
    else
        disp(msg('mpm_incomplete'));
        I = find(~tf).';
        for i = 1:numel(I)
            disp(['- ', reqPackages{I(i)}]);
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
        misPackages = reqPackages(find(~tf).');
        disp('Calling MPM installer -- installing req. packages');
        disp(lineStr);
        for i = 1:numel(misPackages)
            installMissingPackageMPM(misPackages{i}); 
            disp(lineStr);
        end
    end
end
