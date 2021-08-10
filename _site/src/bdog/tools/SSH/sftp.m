function ssh2_struct = sftp(ssh2_struct)
% SFTP   Main SFTP interface for SSH2. The input struct provides all 
%        the parameters for the SSH2 function to perform a sftp command.
%
%   SFTP(SSH2_CONN) enables a ssh2 connection for sftp and runs
%   it through the ssh2 engine where a previous configured scp command
%   will be executed.
% 
%   SFTP returns the SSH2_CONN for future use.
%
%see also sftp_get, sftp_put, sftp_simple_get, sftp_simple_put
%
% (c)2011 Boston University - ECE
%    David Scott Freedman (dfreedma@bu.edu)
%    Version 2.0

if nargin == 1
    ssh2_struct.sftp = 1;
    ssh2_struct = ssh2(ssh2_struct);
else
    ssh2_struct = [];
    help sftp
end
