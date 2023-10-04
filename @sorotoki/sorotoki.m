function vargout = sorotoki(varargin)
%SOROTOKI is the installer function of SOROTOKI.
%   SOROTOKI can be executed in the command window and it will start the
%   installation automatically. Make sure you have an active internet
%   connection to start the installation, which will be compared to the
%   most current (stable) version of the toolkit. Some embedded functions
%   for the SOROTOKI installer are:
%
%   sorotoki install            % installer
%   sorotoki remove             % remove
%   sorotoki disable            % disables toolkit, but keeps source files
%   sorotoki cd                 % go to toolkit
%
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
   
    global log auto_approve
    auto_approve = false;
    % log.options.isDebug = true;

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
    installFile = whereis('sorotoki.m');
    installPath = installFile(1:end-10);

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
        gotoSorotokiFolder(installPath);
        return
    end
    
    if strcmpi(action,'install') || strcmpi(action,'-i')

        if sum(ismember(prompt,'--approve'))
            auto_approve = true;
        end

        installSorotoki(reqPackages,soroPackages);
        addfolder('lib');
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

    % 
    if strcmpi(action,'build') || strcmpi(action,'-b') || ...
       strcmpi(action,'mex')

        if sum(ismember(prompt,'--approve'))
            auto_approve = true;
        end
        buildSorotokiMex(installPath);
        return
    end

    if strcmpi(action,'test') || strcmpi(action,'-t') || ...
        strcmpi(action,'testsuite')
        out = runSorotokiTest(prompt);
        if nargout > 0
            vargout{1} = out;
            assert(vargout{1}, 'One or more tests have failed...');
        end
        return
     end
end