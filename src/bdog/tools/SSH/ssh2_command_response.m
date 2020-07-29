function [command_result] = ssh2_command_response(ssh2_struct)
% SSH2_COMMAND_RESPONSE   Returns the host response from the last SSH2 command
%
%   SSH2_COMMAND(SSH2_CONN)
%   Is a quick method to retrieve the command response from the last SSH2
%   command from the remote host.
% 
%   SSH2_COMMAND returns the COMMAND RESPONSE from the last remote command.
%
%see also ssh2, ssh2_commmand, ssh2_simple_command
%
% (c)2011 Boston University - ECE
%    David Scott Freedman (dfreedma@bu.edu)
%    Version 2.0

if nargout == 1
    if nargin == 1
        command_result = ssh2_struct.command_result;
    else
        command_result = 0;
    end
else
    if (nargout < 2 || enableprint)
        for i = 1:numel(ssh2_struct.command_result)
           fprintf('%s\n', ssh2_struct.command_result{i});
        end
    end
end
