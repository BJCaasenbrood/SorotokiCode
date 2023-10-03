function installMissingPackageMPM(package,varargin)

    action = 'install';
    if ~isempty(varargin)
        action = varargin{1};
    end

    switch package
        case 'matlabessentialskit'
            mpm(action,'matlabessentialskit',...
                '--force','--all-paths','--github-first');
        case 'interparc'
            mpm(action,'interparc','--all-paths');
        case 'inpolyhedron'
            mpm(action,'inpolyhedron','--all-paths');         
        case 'distance2curve'
            mpm(action,'distance2curve','--all-paths');
        case 'matlabgraphicalmodel'             
            mpm(action,'matlabgraphicalmodel','--all-paths','--github-first');
        case 'sorotokisdf'                        
            mpm(action,'sorotokisdf','install-dir','./lib/','--all-paths','--force');
        case 'sorotokimesh'                        
            mpm(action,'sorotokimesh','install-dir','./lib/','--all-paths','--force');     
        case 'sorotokifem'                      
            mpm(action,'sorotokifem','install-dir','./lib/','--all-paths','--force');            
        case 'sorotokimodel'               
            mpm(action,'sorotokimodel','install-dir','./lib/','--all-paths','--force');         
        case 'sorotokibots'            
            if strcmpi(action,'uninstall')
                try
                    movefile assets/* assets/sorotokibots
                end
                mpm(action,'sorotokibots','install-dir','./assets/','--force');  
            else
                mkdir('assets');
                mpm(action,'sorotokibots','install-dir','./assets/','--all-paths','--force','-u','https://github.com/BJCaasenbrood/SorotokiBots.git');   
                try  
                    movefile assets/sorotokibots/* assets
                    addfolder('assets');     
                end
            end
    end
end