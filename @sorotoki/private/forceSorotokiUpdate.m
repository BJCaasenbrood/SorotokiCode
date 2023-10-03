function forceSorotokiUpdate(soroPackages)
    warning off;
    for i = 1:numel(soroPackages)
        installMissingPackageMPM(soroPackages{i}); 
    end
    keyboard;
    warning on;
end