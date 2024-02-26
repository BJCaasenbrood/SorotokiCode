function forceSorotokiUpdate(soroPackages)
    warning off;
    for i = 1:numel(soroPackages)
        installMissingPackageMPI(soroPackages{i}); 
    end
    keyboard;
    warning on;
end