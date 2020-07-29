function ssh2_struct = ssh2_config_publickey(hostname, username, pemFile, pemFilePassword, port)
% SSH2_CONFIG_PUBLICKEY   setups a  SSH2 connection with the
%                         specified hostname, username, and public key
%
%   SSH2_CONFIG_PUBLICKEY(HOSTNAME,USERNAME,PRIVATE_KEY,PRIVATE_KEY_PASSWORD, [PORT])
%   Configures a connection to the host, HOSTNAME with user USERNAME. The
%   private key can be a path to a private key file or a private key string
%   since this is unlikey they can be confused. The password to unencrypt
%   the private key us supplied as PRIVATE_KEY_PASSWORD. The returned
%   SSH2_CONN can be used by SSH2_COMMAND or other SSH2 file transfer
%   commands.
%
%   OPTIONAL INPUTS:
%   -----------------------------------------------------------------------
%   PORT  to specify a non-standard SSH TCP/IP port. Default is 22.
%
% 
%   SSH2_CONFIG_PUBLICKEY returns the SSH2_CONN for future use
%
%see also ssh2_config, ssh2, ssh2_command, scp_get, scp_put, sftp_get, sftp_put
%
% (c)2011 Boston University - ECE
%    David Scott Freedman (dfreedma@bu.edu)
%    Version 2.0

ssh2_struct = ssh2_setup(); %default config

if (nargin >= 4)
    ssh2_struct.hostname = hostname;
    ssh2_struct.username = username;
    if (exist(pemFile,'file'))
        ssh2_struct.pem_file = pemFile;
    else
        ssh2_struct.pem_private_key = pemFile;
    end

    ssh2_struct.pem_private_key_password = pemFilePassword;

    if nargin >= 5
        ssh2_struct.port = port;
    end
else
    help ssh2_config_publickey
end