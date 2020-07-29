function [ssh2_struct, command_result] = ssh2_command(ssh2_struct, command, enableprint)
% SSH2_COMMAND   Reuses a SSH2 connection and issues the specified command
%
%   [SSH2_CONN, [COMMAND_RESULT]] = SSH2_COMMAND(SSH2_CONN,COMMAND,[ENABLEPRINTTOSCREEN])
%   Connects to the SSH2 host with a configured SSH2_CONN. Once connected 
%   the COMMAND is issues. The output from the remote host is returned.
%   The connection to the remote host is still open the function completes.
%   When finished, close the connection with SSH2_CLOSE(SSH2_CONN)
%
%   OPTIONAL INPUTS:
%   -----------------------------------------------------------------------
%   ENABLEPRINTTOSCREEN  force printing the remote host output on screen.
%
% 
%   SSH2_COMMAND returns the SSH2_CONN for future use and the 
%   cell array containing the host response.
%
%see also ssh2, ssh2_config, ssh2_config_publickey, ssh2_simple_command
%
% (c)2011 Boston University - ECE
%    David Scott Freedman (dfreedma@bu.edu)
%    Version 2.0

if nargin < 3
   enableprint = 0;
end

ssh2_struct.command = command;
command_result = [];

%% SSH TO HOST
ssh2_struct = ssh2(ssh2_struct);

%% OUTPUT RESPONSE FROM HOST
if (enableprint)
    for i = 1:numel(ssh2_struct.command_result)
       fprintf('%s \n', ssh2_struct.command_result{i});
    end
end

if nargout > 1
    command_result = ssh2_struct.command_result;
end
if nargout == 0
    clear ssh2_struct;
end
