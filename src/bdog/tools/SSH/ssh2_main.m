function ssh2_struct = ssh2_main(ssh2_struct)
% SSH2_MAIN   main control program for interfacing with Ganymed Java SSH-2 Library
%             SHOULD NOT BE CALLED DIRECTLY, but through ssh2.m
%
%see also ssh2_config, ssh2_config_publickey, ssh2_command, scp_get, scp_put, sftp,get, sftp_put, ssh2_close
%
% (c)2013 Boston University - ECE
%    David Scott Freedman (dfreedma@bu.edu)
%    Version 4.0

%% LOAD JAVA LIB

% http://www.ganymed.ethz.ch/ssh2/javadoc/index.html
%
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.File;

import ch.ethz.ssh2.*;
import ch.ethz.ssh2.Connection;
import ch.ethz.ssh2.Session;
import ch.ethz.ssh2.StreamGobbler;   
import ch.ethz.ssh2.SCPClient;
import ch.ethz.ssh2.SFTPv3Client;
import ch.ethz.ssh2.SFTPv3FileHandle;

%% ASSUME SSH2.M HAS CHECKED SETUP
%

%% INPUTS are OK, PROCEED WITH MAIN FUNCTION
%   will check for lib and connection first, then connect and perform task

% SETUP SSH2 LIBRARY
if (ssh2_struct.ssh2_java_library_loaded == 0)
    error('Error: SSH2 could not load the Ganymed SSH2 java package!'); 
end

%% SETUP CONNECTION

%% AUTO-RECONNECT IF ENABLED
if (ssh2_struct.autoreconnect == 1 && ssh2_struct.authenticated == 1) % assume we're already connected
    try
        fprintf('Checking Connection...\n');
        tmp_command_session  =  ssh2_struct.connection.openSession();
        tmp_command_session.close();
    catch 
        fprintf('Reconnecting...\n');
        ssh2_struct.authenticated = 0;
    end

end

%% AUTHENTICATE (PASSWORD OR PRIVATE KEY)
if (ssh2_struct.authenticated == 1) % assume we're already connected

else
    try %% MAKE A CONNECTION HERE
        ssh2_struct.connection = Connection(ssh2_struct.hostname,...
                                            ssh2_struct.port);
        ssh2_struct.connection.connect();

        if (~isempty(ssh2_struct.pem_private_key))
                           ssh2_struct.authenticated = ssh2_struct.connection.authenticateWithPublicKey(...
            ssh2_struct.username,ssh2_struct.pem_private_key,ssh2_struct.pem_private_key_password);
        elseif (~isempty(ssh2_struct.pem_file))
             try
                 private_key_file_handle = java.io.File(ssh2_struct.pem_file);
             catch 
                 error(['Error: SSH2 could not open private key file'...
                    ' %s ...'],...
                    private_key);
             end
            ssh2_struct.authenticated = ssh2_struct.connection.authenticateWithPublicKey(...
            ssh2_struct.username,private_key_file_handle,ssh2_struct.pem_private_key_password);
        else
            ssh2_struct.authenticated = ssh2_struct.connection.authenticateWithPassword(...
            ssh2_struct.username,ssh2_struct.password);
        end


        if(~ssh2_struct.authenticated)
            error('Error: SSH2 could not authenticate the connection!');
        end
    catch 
        error('Error: SSH2 could not connect to the ssh2 host - "%s"!',...
                ssh2_struct.hostname);
    end
end


%% SCP FILE TRANSFER
if (ssh2_struct.scp > 0) 
    ssh2_struct.scp = 0; % clear for next time
    getFile = 0;
    putFile = 0;
    strArrayGet = [];
    strArrayPut = [];
    strArrayPutRename = [];
    
    scp1 = SCPClient(ssh2_struct.connection);
    if (ssh2_struct.getfiles)
        ssh2_struct.getfiles = 0; %setup for next time
        getFile = 1;
        if (iscell(ssh2_struct.remote_file))
            num_of_files = numel(ssh2_struct.remote_file);
            strArrayGet = javaArray('java.lang.String', num_of_files);
            for remote_file_index = 1:num_of_files
                strArrayGet(remote_file_index) = java.lang.String( ...
                   ssh2_remotePathString(...
                        ssh2_struct.remote_file{remote_file_index},...
                        ssh2_struct.remote_target_direcory) );
            end
        else %assume a string
            strArrayGet = ssh2_remotePathString(...
                                ssh2_struct.remote_file,...
                                ssh2_struct.remote_target_direcory);
        end
        ssh2_struct.remote_file = []; % clear out for next time
    end
    if (ssh2_struct.sendfiles)
        ssh2_struct.sendfiles = 0; %setup for next time
        putFile = 1;
        if (iscell(ssh2_struct.local_file))
            num_of_files = numel(ssh2_struct.local_file);
            strArrayPut = javaArray('java.lang.String', num_of_files);
            for local_file_index = 1:num_of_files
                strArrayPut(local_file_index) = java.lang.String(...
                   ssh2_localPathString(...
                            ssh2_struct.local_file{local_file_index},...
                            ssh2_struct.local_target_direcory) );
            end
        else %assume a string
            strArrayPut = ssh2_localPathString(...
                            ssh2_struct.local_file,...
                            ssh2_struct.local_target_direcory);
        end
        ssh2_struct.local_file = []; % clear out for next time
        if (~isempty(ssh2_struct.remote_file_new_name))
            if (iscell(ssh2_struct.remote_file_new_name))
                num_of_files = numel(ssh2_struct.remote_file_new_name);
                strArrayPutRename = javaArray('java.lang.String', num_of_files);
                for remoteRename_file_index = 1:num_of_files
                    strArrayPutRename(remoteRename_file_index) = java.lang.String(...
                                ssh2_struct.remote_file_new_name{remoteRename_file_index});
                end
            else %assume a string
                strArrayPutRename = ssh2_struct.remote_file_new_name;
            end
        end
        ssh2_struct.remote_file_new_name = []; % clear out for next time
    end
    
    try
        if (getFile > 0)
            scp1.get(strArrayGet,ssh2_struct.local_target_direcory);
        end
        if (putFile > 0)
            if (~isempty(strArrayPutRename)) %rename remote files
                scp1.put(strArrayPut,strArrayPutRename,...
                        ssh2_struct.remote_target_direcory, ...
                        sprintf('%04d',ssh2_struct.remote_file_mode));
            else %use regular file names
                scp1.put(strArrayPut,ssh2_struct.remote_target_direcory, ...
                        sprintf('%04d',ssh2_struct.remote_file_mode));
            end
        end
    catch err
        error('\nError Transferring File\nSEE JAVA ERROR BELOW\n\n%s',err.message);
    end
    clear scp1;
end

%% SFTP FILE TRANSFER
if (ssh2_struct.sftp > 0) 
    ssh2_struct.sftp = 0; % clear for next time
    getFile = 0;
    putFile = 0;
    strArrayGet = [];
    strArrayPut = [];
    strArrayPutRename = [];

    if (ssh2_struct.getfiles)
        ssh2_struct.getfiles = 0; %setup for next time
        getFile = 1;
        if (iscell(ssh2_struct.remote_file))
            strArrayGet = ssh2_struct.remote_file;
        else %assume a string
            strArrayGet{1} = ssh2_struct.remote_file;
        end
        ssh2_struct.remote_file = []; % clear out for next time
    end
    if (ssh2_struct.sendfiles)
        ssh2_struct.sendfiles = 0; %setup for next time
        putFile = 1;
        if (iscell(ssh2_struct.local_file))
            strArrayPut = ssh2_struct.local_file;
        else %assume a string
            strArrayPut{1} = ssh2_struct.local_file;
        end
        ssh2_struct.local_file = []; % clear out for next time
        if (~isempty(ssh2_struct.remote_file_new_name))
            if (iscell(ssh2_struct.remote_file_new_name))
                strArrayPutRename = ssh2_struct.remote_file_new_name;
            else %assume a string
                strArrayPutRename{1} = ssh2_struct.remote_file_new_name;
            end
        end
        ssh2_struct.remote_file_new_name = []; % clear out for next time
    end

    try
        if (getFile > 0)
            %% A custom ganymed-ssh2 (ganymed-ssh2-m1, not ganymed-ssh2-build250)
            %% file enables this functionaly. Its' not efficient to check
            %% whether the custom library is being used because the script will error out
            %% anyways. Instead, a better error message is added to the end.
            
            %% NETWORK TRANSFERS ARE SLOWER, USE SCP IF POSSIBLE
            
            %Open SFTP Connection
            sftp1 = SFTPv3Client(ssh2_struct.connection);
            for getIndex = 1:numel(strArrayGet)
                %create files

                remoteFilePath = ssh2_remotePathString(...
                        strArrayGet{getIndex},...
                        ssh2_struct.remote_target_direcory);
                localFilePath = ssh2_localPathString(...
                        ssh2_remotePathFilenameString(strArrayGet{getIndex}),...
                        ssh2_struct.local_target_direcory);

                remotef=sftp1.openFileRO(remoteFilePath); %read only

                %transfer file byte by byte in 1kB chunks
                count=0;
                readcnt = 0; 
                bufsize = 1024;
                Jbuf = javaArray('java.lang.Byte', bufsize);
                local_fileid = fopen(localFilePath,'w');
                
                
                try
                    while(readcnt~=-1)
                        Jbuf = sftp1.readMatlab(remotef,uint32(count), uint8(zeros(1,bufsize)),uint16(0),uint16(bufsize));
                        readcnt = sftp1.getReadLen();
                        count = count + readcnt;
                        % from Handling Data Returned from Java Methods,
                        % byte arrays will return as signed int8
                        fwrite(local_fileid, Jbuf(1:readcnt),'int8');
                    end   
                catch 
                    extraErrStr = '';
                    if (~strcmp(ssh2_struct.ganymed_java_library,'ganymed-ssh2-m1'))
                        extraErrStr = sprintf('Possible JAVA/MATLAB incompatibility\nCannot use SFTP to retrieve files!\n');
                        extraErrStr = [extraErrStr sprintf(' !! Please ensure the included CUSTOM ganymed-ssh2 library is being used!!\n')];
                        extraErrStr = [extraErrStr sprintf(' see ssh2_setup.m for more information... or use SCP instead!\n')];
                        extraErrStr = [extraErrStr sprintf('See http://www.mathworks.com/matlabcentral/newsreader/view_thread/71084 and\n')];
                        extraErrStr = [extraErrStr sprintf(' http://www.mathworks.com/help/matlab/matlab_external/handling-data-returned-from-a-java-method.html\n for more info.\n\n')];
                    end
                     error(['%sError: SFTP could not write to the file: %s, '...
                        ' on local machine!'],extraErrStr,localFilePath);
                end
                sftp1.closeFile(remotef); %all done!
                fclose(local_fileid); %close the file
            end
            sftp1.close(); %close the connection
        end
        if (putFile > 0)
            %Open SFTP Connection
            sftp1 = SFTPv3Client(ssh2_struct.connection);
            for putIndex = 1:numel(strArrayPut)
                %create files
                if (~isempty(strArrayPutRename)) %rename remote files
                    remoteFilePath = ssh2_remotePathString(...
                        strArrayPutRename{putIndex},...
                        ssh2_struct.remote_target_direcory);
                else
                    remoteFilePath = ssh2_remotePathString(...
                        strArrayPut{putIndex},...
                        ssh2_struct.remote_target_direcory);
                end
                localFilePath = ssh2_localPathString(...
                        strArrayPut{putIndex},...
                        ssh2_struct.local_target_direcory);

                localf=sftp1.createFile(remoteFilePath); %create the file
                remotef=sftp1.openFileRW(remoteFilePath); %get ready to write

               %transfer file byte by byte in 1kB chunks
                count=0; buf = zeros(1,1024);
                local_fileid = fopen(localFilePath,'r');
                [buf, bufsize] = fread(local_fileid, 16*1024);
                try 
                    while(bufsize~=0)
                        sftp1.write(remotef,count,buf,0,bufsize);
                        count=count+bufsize;
                        [buf, bufsize] = fread(local_fileid, 16*1024);
                    end   
                catch
                    error(['Error: SFTP could not write to the file: %s, '...
                        ' on remote machine - "%s"!'],remoteFilePath,ssh2_struct.hostName);
                end
                sftp1.closeFile(remotef); %all done!
                fclose(local_fileid); %close the file
            end
            sftp1.close(); %close the connection
        end
    catch err
        error('\nSFTP: Error Transferring File\nSEE JAVA ERROR BELOW\n\n%s',err.message);
    end
end

%% ISSUE COMMANDS
ssh2_struct.command_result = {''}; % clear this out
if (ischar(ssh2_struct.command)) % we should send a command then.
    % open session and send commands
    ssh2_struct.command_session = ssh2_struct.connection.openSession();
    ssh2_struct.command_session.execCommand(ssh2_struct.command);

    if (~ssh2_struct.command_ignore_response) %get the response from the host
        stdout = StreamGobbler(ssh2_struct.command_session.getStdout());
        br = BufferedReader(InputStreamReader(stdout));
        while(true)
            line = br.readLine();
            if(isempty(line))
                break
            else
                if(isempty(ssh2_struct.command_result{1}))
                    ssh2_struct.command_result{1}  =  char(line);
                else
                    ssh2_struct.command_result{end+1}  =  char(line);
                end
            end
        end
        ssh2_struct.command_result = ssh2_struct.command_result';
    else

    end
    %pause(.1);
    ssh2_struct.command_session.close();
    ssh2_struct.command = []; % clear out the previous command
end

%% CLOSE THE CONNECTION
if(ssh2_struct.close_connection > 0)
    ssh2_struct.connection.close();
    ssh2_struct.authenticated = 0;
end

end

%% CLOSE THE CONNECTION
%if(ssh2_struct.close_connection > 0)
%    ssh2_struct.connection.close();
%    ssh2_struct.authenticated = 0;
%end
    

%% HELPER FUNCTIONS FOR REMOTE DIRECTORY PARSING
function [str] = ssh2_remotePathFilenameString(fileStr)
% for SFTP to get the correct path when downloading from a unix server to
% windows client
[pathstr, name, ext] = fileparts(fileStr);
str = [name ext];

end

function [str] = ssh2_remotePathString(fileStr,pathStr)
% assume that we're uploading to a unix server, filesep = '/';
filesepstr = '/';
if (~isempty(pathStr) && (pathStr(end) ~= filesepstr)) %no pathsep included
    str = [pathStr filesepstr fileStr];
elseif (~isempty(pathStr)) %pathsep is included
    str = [pathStr fileStr];
else
    str = [fileStr];
end

end

function [str] = ssh2_localPathString(fileStr,pathStr)
% use the encoding scheme of local host
filesepstr = filesep();
if (~isempty(pathStr) && (pathStr(end) ~= filesepstr)) %no pathsep included
    str = [pathStr filesepstr fileStr];
elseif (~isempty(pathStr)) %pathsep is included
    str = [pathStr fileStr];
else
    str = [fileStr];
end

end