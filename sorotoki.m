function sorotoki(arg)
%SOROTOKI is the installer function of SOROTOKI.
%   SOROTOKI can be executed in the command window and it will start the
%   installation aumatically. Make sure you have an active internet
%   connection to start the installation, which will be compared to the
%   most current (stable) version of the toolkit. Some embedded functions
%   for the SOROTOKI installer are:
%
%   sorotoki install       % installer
%   sorotoki -i            % ...
%   sorotoki check         % check installation/toolkit
%   sorotoki -c            % check installation/toolkit
%   sorotoki demo          % lists demos
%   sorotoki -d            % lists demos
%   sorotoki help          % help
%   sorotoki -h            % help
%   sorotoki version       % version
%   sorotoki -v            % version
%
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
%
%   Last edit: March 21, 2023.

if nargin < 1, clc; clear; arg = 'install'; end
% ------------------------------------------------------------------------
switch(arg)
    case('install');   installSorotokiToolkit;
    case('uninstall'); unloadSorotoki;
    case('remove');    unloadSorotoki;
    case('demo');      showDemo;

    case('cd');        getPath;
    case('check');     verifySorotoki;
    case('version');   disp(soropatch(1));
    case('--version'); disp(soropatch(1));
    case('-v');        disp(soropatch(1));
    case('-c');        verifySorotoki;
    case('-i');        installSorotokiToolkit;
    case('-h');        help sorotoki;
    case('-d');        showDemo;
    case('-u');        unloadSorotoki;
end
end
% ------------------------------------------------------------------------

function installSorotokiToolkit
addpath('src/__base__/sorotoki/');
installSorotoki;
end
