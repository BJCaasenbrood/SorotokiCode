function removeSorotoki(soroPackages)
    
    warning off;

    disp('Calling MPI installer -- removing Sorotoki');
    for i = 1:numel(soroPackages)
        installMissingPackageMPI(soroPackages{i},'uninstall'); 
    end

    try
        rmdir('lib');
        rmdir('assets');
        delete sorotoki.log
    end

    warning on;
end