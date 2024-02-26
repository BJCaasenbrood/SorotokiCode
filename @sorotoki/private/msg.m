function str = msg(id)
    switch id
        case 'missing_toolboxes'
            str = 'Some prerequisite toolboxes are missing from MATLAB';
        case 'request_to_install_toolbox'        
            str = 'Sorotoki halted! Please install required toolboxes before proceeding';
        case 'mpm_not_installed'
            str = ['MPI is not installed, or overshadowed by build-in mpm! ', ...
                'Please install MPI via the Add-On Manager before proceeding'];
        case 'mpm_package_missing'
            str = 'MPI package not installed';
        case 'mpm_install_missing'
            str = 'Installing missing MPI package from repo';
        case 'mpm_complete'
            str = 'All required package are found on MPI';
        case 'mpm_incomplete'
            str = 'Some required package are missing from MPI';
        case 'soro_complete'
            str = 'Sorotoki libraries found on MPI are up-to-date';
        case 'soro_incomplete'
            str = 'Some Sorotoki packages are missing from MPI';       
    end
end