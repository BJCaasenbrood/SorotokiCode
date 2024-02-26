function setSorotokiSource(soroPackages,opt)
    warning off;
    for i = 1:numel(soroPackages)-1
        installMissingPackageMPI(soroPackages{i},opt); 
    end
    warning on;
end