function ssh2_struct = sftp_get(ssh2_struct, remoteFilename, localPath, remotePath)
%  SFTP_GET   creates a simple SSH2 connection and retrieves a remote file
%             NOTE: BECAUSE OF THE WAY MATLAB WORKS WITH JAVA BYTE ARRAYS
%             GARYMED-SSH2 is CUMBERSOME TO USE WITH SFTP.
%             SCP IS USED (INPLACE) INSTEAD
% SFTP_GET   Reuse configured ssh2_connection to SFTP remote files to local host
%
%   SFTP_GET(SSH2_CONN,REMOTEFILENAME,[LOCALPATH],[REMOTEPATH])
%   uses a ssh2_connection and downloads the REMOTEFILENAME from the remote
%   host using SFTP. SSH2_CONN must already be confgured using 
%   ssh2_config or ssh2_config_publickey.
%
%   REMOTEFILENAME can be either a single string, or a cell array of strings. 
%   If REMOTEFILENAME is a cell array, all files will be downloaded
%   sequentially.
%
%   OPTIONAL INPUTS:
%   -----------------------------------------------------------------------
%   LOCALPATH specifies the folder to download the remote file to.
%   Otherwise the working directory is used.
%   REMOTEPATH specifies a specific path on the remote host to look for the 
%   file to download. Otherwise, the default (home) folder is used.
%
%   SFTP_GET returns the SSH2_CONN for future use.
%
%see also sftp_put, sftp_simple_get, sftp_simple_put, sftp
%
% (c)2013 Boston University - ECE
%    David Scott Freedman (dfreedma@bu.edu)
%    Version 4.0

if nargin < 2
    if nargin == 0
        ssh2_struct = [];
    end
    help sftp_get
else
    if nargin < 3
        localPath = pwd();
    elseif isempty(localPath)
        localPath = pwd();   
    end    
    
    if nargin < 4
        remotePath = '';         
    end
    
    ssh2_struct.getfiles = 1;
    ssh2_struct.remote_file = remoteFilename;
    ssh2_struct.local_target_direcory = localPath;
    ssh2_struct.remote_target_direcory = remotePath;
    %ssh2_struct.local_file = localFilename; unused in SFTP, filename will
    %be taken from remoteFilename

    ssh2_struct = sftp(ssh2_struct);
end

