function ssh2_struct = ssh2_config(hostname, username, password, port)
% SSH2_CONFIG   creates a simple SSH2 connection with the
%               specified hostname, username, and password
%
%   SSH2_CONFIG(HOSTNAME,USERNAME,PASSWORD, [PORT])
%   Configures a connection to the host, HOSTNAME with user USERNAME and
%   password, PASSWORD. The returned SSH2_CONN can be used by SSH2_COMMAND 
%   or other SSH2 file transfer commands.
%
%   OPTIONAL INPUTS:
%   -----------------------------------------------------------------------
%   PORT  to specify a non-standard SSH TCP/IP port. Default is 22.
%
% 
%   SSH2_CONFIG returns the SSH2_CONN for future use
%
%see also ssh2_config_publickey, ssh2, ssh2_command, scp_get, scp_put, sftp_get, sftp_put
%
% (c)2011 Boston University - ECE
%    David Scott Freedman (dfreedma@bu.edu)
%    Version 2.0

ssh2_struct = ssh2_setup(); %default config
if (nargin >= 3)
    ssh2_struct.hostname = hostname;
    ssh2_struct.username = username;
    ssh2_struct.password = password;
    if nargin >= 4
        ssh2_struct.port = port;
    end
else
    help ssh2_config
end