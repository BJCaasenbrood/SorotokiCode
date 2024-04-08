function checkSoroPackages(soroPackages,mpiPackages)

    global auto_approve

    disp('Checking for updates: Sorotoki library'); 
    installed = ismember(soroPackages,mpiPackages);
    
    if all(installed)

        % retrieve date MPI information
        [~, mpiDates] = getPackagesInformationMPI();

        % check if packages are outdated
        outdated = [];
        for i = 1:numel(soroPackages)
            [~, lastestDate] = getCommitInformationGIT(soroPackages{i})
            I = find(ismember(mpiPackages,soroPackages{i}));

            installDate = datetime(mpiDates{I}, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss')
            if installDate >= lastestDate
                outdated = [outdated; false];
            else
                outdated = [outdated; true];
            end
        end

        if all(~outdated)
            disp(msg('soro_complete'));
        else
            disp(msg('soro_outofdate'));
            I = find(~outdated).';
            for i = 1:numel(I)
                disp(['- ', soroPackages{I(i)}]);
            end

            if ~auto_approve
                disp('Update MPM packages?');
                reply = input(i18n('confirm'), 's');
                if isempty(reply)
                    reply = i18n('confirm_yes');
                end
                if ~strcmpi(reply(1), i18n('confirm_yes'))
                    disp(i18n('confirm_nvm'));
                    return;
                end
            end
        
            lineStr = repmat('━', 1, 40);
            misPackages = soroPackages(find(~outdated).');
            disp('Calling MPI installer -- updating Sorotoki packages');
            disp(lineStr);
            for i = 1:numel(misPackages)
                installMissingPackageMPI(misPackages{i}); 
                disp(lineStr);
            end
        end
    else
        disp(msg('soro_incomplete'));
        I = find(~installed).';
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
                disp(i18n('confirm_nvm'));
                return;
            end
        end
    
        lineStr = repmat('━', 1, 40);
        misPackages = soroPackages(find(~installed).');
        disp('Calling MPI installer -- installing Sorotoki packages');
        disp(lineStr);
        for i = 1:numel(misPackages)
            installMissingPackageMPI(misPackages{i}); 
            disp(lineStr);
        end
    end
end
