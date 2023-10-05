function removeSorotoki(soroPackages)
    
    warning off;

    disp('Calling MPM installer -- removing Sorotoki');
    for i = 1:numel(soroPackages)
        installMissingPackageMPM(soroPackages{i},'uninstall'); 
    end

    try
        rmdir('lib');
        rmdir('assets');
        delete sorotoki.log
    end

    warning on;
end