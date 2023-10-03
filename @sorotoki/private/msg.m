function str = msg(id)
    switch id
        case 'missing_toolboxes'
            str = 'Some prerequisite toolboxes are missing from MATLAB';
        case 'request_to_install_toolbox'        
            str = 'Sorotoki halted! Please install required toolboxes before proceeding';
        case 'mpm_not_installed'
            str = ['MPM is not installed, or overshadowed by build-in mpm! ', ...
                'Please install MPM via the Add-On Manager before proceeding'];
        case 'mpm_package_missing'
            str = 'MPM package not installed';
        case 'mpm_install_missing'
            str = 'Installing missing MPM package from repo';
        case 'mpm_complete'
            str = 'All required package are found on MPM';
        case 'mpm_incomplete'
            str = 'Some required package are missing from MPM';
        case 'soro_complete'
            str = 'Sorotoki libraries found on MPM are up-to-date';
        case 'soro_incomplete'
            str = 'Some Sorotoki packages are missing from MPM';       
    end
end