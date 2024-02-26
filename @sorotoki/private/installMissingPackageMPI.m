function installMissingPackageMPI(package,varargin)

    action = 'install';
    if ~isempty(varargin)
        action = varargin{1};
    end

    switch package
        case 'matlabessentialskit'
            mpi(action,'matlabessentialskit',...
                '--force','--all-paths','--github-first');
        case 'interparc'
            mpi(action,'interparc','--all-paths');
        case 'inpolyhedron'
            mpi(action,'inpolyhedron','--all-paths');         
        case 'distance2curve'
            mpi(action,'distance2curve','--all-paths');
        case 'matlabgraphicalmodel'             
            mpi(action,'matlabgraphicalmodel','--all-paths','--github-first');
        case 'sorotokisdf'                        
            mpi(action,'sorotokisdf','install-dir','./lib/','--all-paths','--force');
        case 'sorotokimesh'                        
            mpi(action,'sorotokimesh','install-dir','./lib/','--all-paths','--force');     
        case 'sorotokifem'                      
            mpi(action,'sorotokifem','install-dir','./lib/','--all-paths','--force');            
        case 'sorotokimodel'               
            mpi(action,'sorotokimodel','install-dir','./lib/','--all-paths','--force');       
        case 'sorotokicontrol'               
            mpi(action,'sorotokicontrol','-u','https://github.com/BJCaasenbrood/SorotokiControl.git',...
                'install-dir','./lib/','--all-paths','--force');                   
        case 'sorotokibots'            
            if strcmpi(action,'uninstall')
                try
                    movefile assets/* assets/sorotokibots
                end
                mpi(action,'sorotokibots','install-dir','./assets/','--force');  
            else
                mkdir('assets');
                mpi(action,'sorotokibots','install-dir','./assets/','--all-paths','--force','-u','https://github.com/BJCaasenbrood/SorotokiBots.git');   
                try  
                    movefile assets/sorotokibots/* assets
                    addfolder('assets');     
                end
            end
    end
end