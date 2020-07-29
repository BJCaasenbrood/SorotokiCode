function ssh2_struct = ssh2_close(ssh2_struct)
% SSH2_CLOSE   closes an existing SSH2 connection. 
%              This function can be used as the last SSH2 commmand
%              because the connection won't be closed until the very end.
%   SSH2_CLOSE(SSH2_CONN) takes a ssh2_connection, sets it close, and runs
%   it through the ssh2 engine where it is closed at the end. 
% 
%   SSH2_CLOSE returns the SSH2_CONN for future use.
%
%see also ssh2, ssh2_setup, ssh2_config, ssh2_config_publickey, scp, sftp
%
% (c)2011 Boston University - ECE
%    David Scott Freedman (dfreedma@bu.edu)
%    Version 2.0

if nargin == 1
    ssh2_struct.close_connection = 1;%;
    ssh2_struct = ssh2(ssh2_struct); %close the connection
else
    ssh2_struct = [];
    help ssh2_close
end