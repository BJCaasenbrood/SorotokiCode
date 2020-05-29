function varargout = git(varargin)
% GIT Execute a git command.
%
% GIT <ARGS>, when executed in command style, executes the git command and
% displays the git outputs at the MATLAB console.
%
% STATUS = GIT(ARG1, ARG2,...), when executed in functional style, executes
% the git command and returns the output status STATUS.
%
% [STATUS, CMDOUT] = GIT(ARG1, ARG2,...), when executed in functional
% style, executes the git command and returns the output status STATUS and
% the git output CMDOUT.

% Check output arguments.
nargoutchk(0,2)

% Specify the location of the git executable.
gitexepath = "C:\Program Files\Git\bin\git.exe";

% Construct the git command.
cmdstr = strjoin([gitexepath, varargin]);

% Execute the git command.
[status, cmdout] = system(cmdstr);
    
switch nargout
    case 0
        disp(cmdout)
    case 1
        varargout{1} = status;
    case 2
        varargout{1} = status;
        varargout{2} = cmdout;
end

end