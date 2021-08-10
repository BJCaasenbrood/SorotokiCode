function command_result = ssh2_simple_command(hostname, username, password, command, enableprint)
% SSH2_SIMPLE_COMMAND   connects to host using given username and password.
%                       A command can be given and the output will be
%                       displayed unless supressed.
%
%   SSH2_SIMPLE_COMMAND(HOSTNAME,USERNAME,PASSWORD,COMMAND,[ENABLEPRINTTOSCREEN])
%   Connects to the SSH2 host, HOSTNAME with supplied USERNAME and
%   PASSWORD. Once connected the COMMAND is issues. The output from the 
%   remote host is returned and the connection is closed.
%
%   OPTIONAL INPUTS:
%   -----------------------------------------------------------------------
%   ENABLEPRINTTOSCREEN  force printing the remote host output on screen.
%
% 
%   SSH2_SIMPLE_COMMAND returns a cell array containing the host response.
%
%see also ssh2_command, ssh2, scp_simple_get, scp_simple_put, sftp_simple_get, sftp_simple_put
%
% (c)2011 Boston University - ECE
%    David Scott Freedman (dfreedma@bu.edu)
%    Version 2.0

command_result = [];

if nargin < 4
    ssh2_struct = [];
    help ssh2_simple_command
else
    if nargin < 5
        enableprint = 0;
    end

    ssh2_struct = ssh2_config(hostname, username, password);
    ssh2_struct.close_connection = 1; %close connection use

    if nargout == 0
        ssh2_struct = ssh2_command(ssh2_struct, command, enableprint);
    else
        [ssh2_struct, command_result] = ssh2_command(ssh2_struct, command, enableprint);
    end
end
