function ssh2_struct = sftp_simple_get(hostname, username, password, remoteFilename, localPath, remotePath)
% SFTP_SIMPLE_GET   creates a simple SSH2 connection and retrieves a remote file
%                   NOTE: BECAUSE OF THE WAY MATLAB WORKS WITH JAVA BYTE ARRAYS
%                   GARYMED-SSH2 is CUMBERSOME TO USE WITH SFTP GETS.
%                   SCP IS USED INSTEAD
%
%   SFTP_SIMPLE_GET(HOSTNAME,USERNAME,PASSWORD,REMOTEFILENAME,[LOCALPATH],[REMOTEPATH])
%   Connects to the SSH2 host, HOSTNAME with supplied USERNAME and
%   PASSWORD. Once connected the REMOTEFILENAME is downloaded from the
%   remote host using SFTP. The connection is then closed.
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
%   SFTP_SIMPLE_GET returns the SSH2_CONN for future use.
%
%see also sftp_simple_put, sftp, sftp_get, sftp_put, ssh2_simple_command
%
% (c)2011 Boston University - ECE
%    David Scott Freedman (dfreedma@bu.edu)
%    Version 2.0

if nargin < 4
    ssh2_struct = [];
    help sftp_simple_get
else
    if nargin < 5
        localPath = pwd();
    elseif isempty(localPath)
        localPath = pwd();   
    end    
    
    if nargin < 6
        remotePath = '';         
    end
    

    ssh2_struct = ssh2_config(hostname, username, password);
    ssh2_struct.close_connection = 1; %close connection use
    ssh2_struct = sftp_get(ssh2_struct, remoteFilename, localPath, remotePath);
end

