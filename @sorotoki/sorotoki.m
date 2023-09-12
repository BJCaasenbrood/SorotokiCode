function sorotoki(varargin)
%SOROTOKI is the installer function of SOROTOKI.
%   SOROTOKI can be executed in the command window and it will start the
%   installation aumatically. Make sure you have an active internet
%   connection to start the installation, which will be compared to the
%   most current (stable) version of the toolkit. Some embedded functions
%   for the SOROTOKI installer are:
%
%   sorotoki install            % installer
%   sorotoki remove             % remove
%   sorotoki disable            % disables toolkit, but keeps source files
%   sorotoki cd                 % go to toolkit
%%
%   Also, some additional toolboxes are required:
%
%   - Image Processing Toolbox
%   - Partial Differential Equation Toolbox
%   - MATLAB Coder
%   - Image Processing Toolbox
%
%   The installer will automaticall check for these. Furthermore, it will
%   also prompt for configurating a startup file. This will ensure that
%   sorotoki is also loaded when opening MATLAB, and can be found under
%   .../MATLAB/startup.m in the MATLAB's main folder on your desktop.
   
    soroPackages = {};
    reqToolboxes = {};
    reqPackages  = {};

    % reading list of requirements for installation
    FID = fopen('requirements.txt','r');
    rdl = textscan(FID,'%q','Delimiter','\n');
    req = rdl{1};
    fclose(FID);

    for ii = 1:numel(req)
        if ~isempty(strfind(req{ii}, '#'))
            % skip commented 
        elseif ~isempty(strfind(req{ii}, 'sorotoki'))
            soroPackages{end+1} = req{ii};
        elseif ~isempty(strfind(req{ii}, ' '))
            reqToolboxes{end+1} = req{ii};
        else
            reqPackages{end+1} = req{ii};
        end
    end

    action = '-h';
    prompt = [];
    mpmPath = which('mpm');
    assert(~strcmpi(mpmPath(1:5),'build'), msg('mpm_not_installed'));
    checkToolboxes(reqToolboxes);
    
    if ~isempty(varargin)
        action = varargin{1};
        if numel(varargin) > 1
            prompt = {varargin{2:end}};
        end
    end

    if strcmpi(action,'help') || strcmpi(action,'-h')
        help sorotoki
        return
    end

    if strcmpi(action,'cd') 
        installPath = whereis('sorotoki.log');
        %installPath = logFile(1:end-12);
        cd(installPath);
        return
    end
    
    if strcmpi(action,'install') || strcmpi(action,'-i')
         installSorotoki(reqPackages,soroPackages);
         addfolder('src');
         return
    end

    if strcmpi(action,'remove') || strcmpi(action,'-r')
        if isempty(prompt)
            removeSorotoki(soroPackages);
        else
            removeSorotoki(prompt);
        end
        return
   end

    if strcmpi(action,'update') || strcmpi(action,'-u')
        if isempty(prompt)
            forceSorotokiUpdate(soroPackages);
        else
            forceSorotokiUpdate(prompt);
        end
        return
    end

    if strcmpi(action,'disable') || strcmpi(action,'-u')
        setSorotokiSource(soroPackages,'disable');
        return
    end

    if strcmpi(action,'build') || strcmpi(action,'-b') || ...
       strcmpi(action,'mex')
        installPath = whereis('sorotoki.log');
        buildSorotokiMex(installPath);
        return
    end

end

% -------------------------------------------------------------------------
function installSorotoki(reqPackages,soroPackages)
    
    installDir = pwd;
    mpmPath    = which('mpm');
    [mpmPackages, mpmPackagedate] = getPackagesInformationMPM();

    if ~isfile('sorotoki.log')
        delete sorotoki.log

        diary sorotoki.log
        disp(['Installation dir: ', installDir]);
        disp(['Directory MPM found: ', mpmPath(1:end-6)]);
    end
    
    checkPackagesMPM(reqPackages,mpmPackages);
    checkSoroPackages(soroPackages,mpmPackages);
    
    diary off;
end
% -------------------------------------------------------------------------
function forceSorotokiUpdate(soroPackages)
    warning off;
    for i = 1:numel(soroPackages)
        installMissingPackgeMPM(soroPackages{i}); 
    end
    warning on;
end
% -------------------------------------------------------------------------
function setSorotokiSource(soroPackages,opt)
    warning off;
    for i = 1:numel(soroPackages)-1
        installMissingPackgeMPM(soroPackages{i},opt); 
    end
    warning on;
end
% -------------------------------------------------------------------------
function removeSorotoki(soroPackages)
    warning off;
    for i = 1:numel(soroPackages)
        installMissingPackgeMPM(soroPackages{i},'uninstall'); 
    end

    try
        rmdir('lib');
        rmdir('assets');
        delete sorotoki.log
    end
    warning on;
end
% -------------------------------------------------------------------------
function buildSorotokiMex(installPath)
    path = installPath;

    disp('Build Sorotoki .mex executables?');
    reply = input(i18n('confirm'), 's');
    if isempty(reply)
        reply = i18n('confirm_yes');
    end
    if ~strcmpi(reply(1), i18n('confirm_yes'))
        disp(i18n('confirm_nvm'));
        return;
    end

end

% -------------------------------------------------------------------------
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
        str = 'Some required package from Sorotoki are missing from MPM';       
end
end

% -------------------------------------------------------------------------
function checkToolboxes(reqToolboxes)

    v = ver;
    installedToolboxes = {v.Name};
    missingToolboxes   = {};
    
    for i = 1:numel(reqToolboxes)
        if ~ismember(reqToolboxes{i},installedToolboxes)
            if isempty(missingToolboxes)
                disp(msg('missing_toolboxes'));
            end
            disp([' - ', reqToolboxes{i}]);
            missingToolboxes{end+1} = reqToolboxes{i};
        end
    end
    
    if ~isempty(missingToolboxes)
        assert(false, msg('request_to_install_toolbox'))
    end

end
% -------------------------------------------------------------------------
function checkPackagesMPM(reqPackages,mpmPackages)

    disp('Checking for updates: MPM package library'); 
    tf = ismember(reqPackages,mpmPackages);

    if all(tf)
        disp(msg('mpm_complete'));
    else
        disp(msg('mpm_incomplete'));
        for i = find(~tf).'
            disp([' - ', reqPackages{i}]);
        end
    
        disp('Install missing MPM packages?');
        reply = input(i18n('confirm'), 's');
        if isempty(reply)
            reply = i18n('confirm_yes');
        end
        if ~strcmpi(reply(1), i18n('confirm_yes'))
            disp(i18n('confirm_nvm'));
            return;
        end
    
        misPackages = reqPackages(find(~tf).');
        for i = 1:numel(misPackages)
            installMissingPackgeMPM(misPackages{i}); 
        end
    end
end

function checkSoroPackages(soroPackages,mpmPackages)

    disp('Checking for updates: Sorotoki library'); 
    tf = ismember(soroPackages,mpmPackages);
    
    if all(tf)
        disp(msg('soro_complete'));
    else
        disp(msg('soro_incomplete'));
        for i = find(~tf).'
            disp([' - ', soroPackages{i}]);
        end
    
        disp('Install missing MPM packages?');
        reply = input(i18n('confirm'), 's');
        if isempty(reply)
            reply = i18n('confirm_yes');
        end
        if ~strcmpi(reply(1), i18n('confirm_yes'))
            disp(i18n('confirm_nvm'));
            return;
        end
    
        misPackages = soroPackages(find(~tf).');
        for i = 1:numel(misPackages)
            installMissingPackgeMPM(misPackages{i}); 
        end
    end
end

% -------------------------------------------------------------------------
function installMissingPackgeMPM(package,varargin)

    action = 'install';
    if ~isempty(varargin)
        action = varargin{1};
    end

    switch package
        case 'matlabessentialskit'
            mpm(action,'matlabessentialskit',...
                '--force','--all-paths','--github-first');
        case 'matlabprogressbar'
            mpm(action,'matlabprogressbar','--all-paths','--force');
        case 'interparc'
            mpm(action,'interparc','--all-paths');
        case 'inpolyhedron'
            mpm(action,'inpolyhedron','--all-paths');         
        case 'distance2curve'
            mpm(action,'distance2curve','--all-paths');
        case 'matlabgraphicalmodel'             
            mpm(action,'matlabgraphicalmodel','--all-paths');
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
% -------------------------------------------------------------------------
function [name, date] = getPackagesInformationMPM()
    mpmPath = which('mpm');
    mpmPackageFolder = [mpmPath(1:end-6),'/mpm-packages/'];
    MPM = load([mpmPackageFolder,'mpm.mat']);
    name = {MPM.packages.name};
    date = {MPM.packages.downloadDate};
end

% -------------------------------------------------------------------------
function str = i18n(key, varargin)
    persistent locale nls
    if ~isstruct(nls)
        locale = char(regexp(get(0, 'Language'), '^[a-zA-Z]+', 'match'));
        nls = load('mpm_nls.mat');
    end

    %% Check if message key exists.
    if isfield(nls, locale)
        data = nls.(locale);
    else
        data = nls.en;
    end

    if ~ischar(key)
        if (                                                                ...
            ~isfield(nls, locale)                                           ...
            || ~isfield(nls.(locale), 'unexpected_key')                     ...
        )
            error(nls.en.unexpected_key, class(key));
        else
            error(nls.(locale).unexpected_key, class(key));
        end
    end

    if (                                                                    ...
        ~isfield(data, key)                                                 ...
        && strcmp(locale, 'en')                                             ...
        || ~isfield(nls.en, key)                                            ...
    )
        error(nls.en.undefined_key, key);
    end

    %% Get the localised message.
    if isfield(data, key)
        str = data.(key);
    else
        str = sprintf(nls.en.(key), string(varargin{:}));
    end

    %% Variable argument substitution.
    if nargin == 1
        return;
    end

    xs = cellfun(@char, varargin, 'uni', false);
    str = sprintf(str, xs{:});
end
