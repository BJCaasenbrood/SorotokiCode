function ssh2_struct = ssh2_setup(ssh2_struct)
% SSH2_SETUP   Create configuration structure, SSH2_CONN,
%               for use with ssh2, scp, sftp and supporting functions. 
%               The file must be created, either directly or indirectly,
%               in order to create a ssh connection. 
%
%               SSH2_CONFIG or SSH2_SIMPLE commands will 
%               automatically create this file for you.
%
% SSH2_SETUP([SSH2_CONN])  normally, no inputs are provided.
%
% In addition to creating (or checking) a default SSH2_CONN config to use, 
% this file will also perform the following steps to ensure the ganymed-ssh2
% java library is loaded in MATLAB.
% 
% If you would like to use SFTP-GET, it is recommended that you change 
% the variable USE_CUSTOM_GANYMED_LIB = 1 %% (line 107)
% In this case, an included custom ganymed-ssh2 library
% will be used that supports this function. The default value of 
% USE_CUSTOM_GANYMED_LIB is 0 in which 1) is performed below
% -----------------------------------------------------------------------
%
% 1.  Download and unzip the ganymed-ssh2 jar file from (if necessary)
%                       http://www.cleondris.ch/opensource/ssh2/
% 2.  Add the jar file to the Matlab dynamic or static Java path
%                       e.g.  to add as a dynamic library 
%                             issue the command
%
%                             javaaddpath('ganymed-ssh2-build250.jar')
%
%                             assuming the SSH-2 Java library is in your
%                             working directory and named
%                             ganymed-ssh2-build250.jar
%
%   OPTIONAL INPUTS:
%   -----------------------------------------------------------------------
%   SSH2_CONN  when an input is supplied, the SSH2_CONN struct will be 
%              examined to verify all necessary fields are present.
%
%
% SSH2_CONN Fields:
% ------------------------------------------------------------------------
%
%   CONNECTION VALUES
%   .hostname  - network name of the SSH2 host
%   .username  - username to login to host as
%   .password  - password associated with username on remote host
%   .port  - SSH2 remote TCP/IP port, default is 22
%
%   PUBLIC/PRIVATE KEY VALUES
%   .pem_file  - file location of the private key
%   .pem_private_key  - private key as a string
%   .pem_private_key_password  - password for the private key
%
%   REMOTE COMMAND PARAMETERS AND OPTIONS
%   .command  - command line string to issue to remote host
%   .command_session  - session object for the ssh2 command line 
%   .command_ignore_response  - set to 1 to ignore response from host
%   .command_result  - cell array containing the response from the host
%
%   GANYMED JAVA OBJECT
%   .connection  - java connection object
%   .ssh2_java_library_loaded  - flag to determine whether ganymed lib is loaded
%
%   CONNECTION FLAGS AND OPTIONS
%   .authenticated  - flag to determine whether we've authenticated
%   .autoreconnect  - flag to check whether or not to reconnect if disconnected
%   .close_connection  - flag to close connection at the end of ssh2.m
%
%   SCP and SFTP FLAGS
%   .scp  - set to 1 to enable an SCP transfer
%   .sftp  - set to 1 to enable an SFTP transfer
%   .sendfiles  - set to 1 to send a file to the remote host
%   .getfiles  - set to 1 to download a file from the remote host
%
%   FILE AND PATHS, LOCAL AND REMOTE
%   .remote_file  - string or cell array of remote files to up/download
%   .local_target_direcory  - string of path to local files to up/download
%   .local_file  - string or cell array of local files to up/download
%   .remote_target_direcory  - string of path to remote files to up/download
%   .remote_file_new_name  - string or cell array of new names of remote files
%   .remote_file_mode  - Integer specifying new uploaded files permissions, default is 0600
%
%   GANYMED VALUES
%   .ganymed_java_library  - name of the current ganymed java library
%   .ganymed_java_library_zip  - filename of ganymed lib to download
%   .ganymed_java_library_jar  - filename of ganymed java lib
%   .ganymed_java_library_jar_path  - path to ganymed Jar file
%
%   UNUSED
%   .verified_config  - Unused variable to bypass config check in ssh2.m
%
%         
%   SSH2_SETUP returns the SSH2_CONN for future use.
%
%see also ssh2_config, ssh2_config_publickey, ssh2, ssh2_simple_command
%
% (c)2013 Boston University - ECE
%    David Scott Freedman (dfreedma@bu.edu)
%    Version 4.0


%% SETUP SSH2 CONNECTION STRUCT
% now includes custom ganymed library, auto download code is still default

% 0 for standard lib, 1 for custom (SFTP-GET support) ganymed-ssh2 library
USE_CUSTOM_GANYMED_LIB = 1; 
% CHANGING THE VARIABLE USE_CUSTOM_GANYMED_LIB REQUIRES RESTARTING MATLAB OR 
% CLEARING previous java DYNAMIC PATH if SSH2 has already been run this
% session.
% NETWORK TRANSFERS ARE SLOWER, USE SCP IF POSSIBLE


if (USE_CUSTOM_GANYMED_LIB > 0)
    ganymed_java_library = 'ganymed-ssh2-m1'; %included custom ganymed library
     ganymed_java_library_zip = [ganymed_java_library '.zip'];
    ganymed_java_library_http = '';
else
    %% 2018-08-20, this has been deprecaited as the library is no longer available from a central website
    %%             instead it is included with this zip file
    %%
    % PICK EITHER CLEONDRIS or GOOGLECODE, NOT BOTH
    % FROM CLEONDRIS
    ganymed_java_library = 'ganymed-ssh2-build250'; % official download from cleondris.ch
    ganymed_java_library_zip = [ganymed_java_library '.zip'];
    ganymed_java_library_http = sprintf('http://www.cleondris.ch/ssh2/%s',ganymed_java_library_zip);
    
    % FROM GOOGLECODE
    %ganymed_java_library = 'ganymed-ssh2-build251beta1'; % official download from cleondris.ch
    %ganymed_java_library_zip = [ganymed_java_library '.zip'];
    %ganymed_java_library_http = sprintf('http://ganymed-ssh-2.googlecode.com/files/%s',ganymed_java_library_zip);
end


ganymed_java_library_jar = [ganymed_java_library '.jar'];
ganymed_java_library_jar_path = [ganymed_java_library ...
                                    filesep() ganymed_java_library_jar];
                                
ssh_path = fileparts(mfilename('fullpath')); 
ganymed_java_library = fullfile(ssh_path,ganymed_java_library); 
ganymed_java_library_zip = fullfile(ssh_path,ganymed_java_library_zip); 
ganymed_java_library_jar = fullfile(ssh_path,ganymed_java_library_jar); 
ganymed_java_library_jar_path = fullfile(ssh_path,ganymed_java_library_jar_path);                                
                                
error_message = 0;
if nargin == 0 %SETUP THE DEFAULT CONFIG
    ssh2_struct.ganymed_java_library = ganymed_java_library;
    ssh2_struct.ganymed_java_library_zip = ganymed_java_library_zip;
    ssh2_struct.ganymed_java_library_jar = ganymed_java_library_jar;
    ssh2_struct.ganymed_java_library_jar_path = ganymed_java_library_jar_path;
    
    ssh2_struct.hostname = [];
    ssh2_struct.username = [];
    ssh2_struct.password = [];
    ssh2_struct.port = 22;

    ssh2_struct.connection = [];
    ssh2_struct.authenticated = 0;
    ssh2_struct.autoreconnect = 0;
    ssh2_struct.close_connection = 0;
    
    ssh2_struct.pem_file = [];
    ssh2_struct.pem_private_key = [];
    ssh2_struct.pem_private_key_password = [];
    
    ssh2_struct.command = [];
    ssh2_struct.command_session = [];
    ssh2_struct.command_ignore_response = 0;
    ssh2_struct.command_result = [];
    
    ssh2_struct.sftp = 0;
    ssh2_struct.scp = 0;
    ssh2_struct.sendfiles = 0;
    ssh2_struct.getfiles = 0;
    
    ssh2_struct.remote_file = [];
    ssh2_struct.local_target_direcory = [];
    ssh2_struct.local_file = [];
    ssh2_struct.remote_target_direcory = [];
    ssh2_struct.remote_file_new_name = [];
    ssh2_struct.remote_file_mode = 0600; %0600 is default
    
    ssh2_struct.verified_config = 0;
    ssh2_struct.ssh2_java_library_loaded = 0;
else
    error_message = 1;
    if (isstruct(ssh2_struct))
        if (isfield(ssh2_struct,'ganymed_java_library') && ...
            isfield(ssh2_struct,'ganymed_java_library_zip') && ...
            isfield(ssh2_struct,'ganymed_java_library_jar') && ...
            isfield(ssh2_struct,'ganymed_java_library_jar_path') && ...
            isfield(ssh2_struct,'hostname') && ...
            isfield(ssh2_struct,'username') && ...
            isfield(ssh2_struct,'password') && ...
            isfield(ssh2_struct,'port') && ...
            isfield(ssh2_struct,'connection') && ...
            isfield(ssh2_struct,'authenticated') && ...
            isfield(ssh2_struct,'autoreconnect') && ...
            isfield(ssh2_struct,'close_connection') && ...
            isfield(ssh2_struct,'pem_file') && ...
            isfield(ssh2_struct,'pem_private_key') && ...
            isfield(ssh2_struct,'pem_private_key_password') && ...
            isfield(ssh2_struct,'command') && ...
            isfield(ssh2_struct,'command_session') && ...
            isfield(ssh2_struct,'command_ignore_response') && ...
            isfield(ssh2_struct,'command_result') && ...
            isfield(ssh2_struct,'sftp') && ...
            isfield(ssh2_struct,'scp') && ...
            isfield(ssh2_struct,'sendfiles') && ...
            isfield(ssh2_struct,'getfiles') && ...
            isfield(ssh2_struct,'remote_file') && ...
            isfield(ssh2_struct,'local_target_direcory') && ...
            isfield(ssh2_struct,'local_file') && ...
            isfield(ssh2_struct,'remote_target_direcory') && ...
            isfield(ssh2_struct,'remote_file_new_name') && ...
            isfield(ssh2_struct,'remote_file_mode') && ...
            isfield(ssh2_struct,'verified_config') && ...
            isfield(ssh2_struct,'ssh2_java_library_loaded') )
            error_message = 0;
        end
    end
    if (error_message == 1)
        warning('SSH2_SETUP: Invalid input provided!');
        help ssh2_setup
    end        
end

%% LOAD GANYMED LIBRARY IF NECESSARY
%   will automatically download the library if it's not available
if (error_message == 0 && ~ssh2_struct.ssh2_java_library_loaded)
    [ssh2_struct] = ssh2_config_check_for_java_lib(ssh2_struct);
    if (~ssh2_struct.ssh2_java_library_loaded)
        if (exist(ganymed_java_library_jar_path,'file'))
            fprintf('Adding Ganymed-ssh2 to the java path by running\njavaaddpath(''%s'')\n',ganymed_java_library_jar_path);
            javaaddpath(ganymed_java_library_jar_path);
            fprintf('\nJust added Ganymed-ssh2 to Matlab''s dynamic java Classpath.\n');
        else
            if (USE_CUSTOM_GANYMED_LIB == 0) %try to download official (non-custom) version 
                if (exist(ganymed_java_library_zip))
                    % FILE IS INCLUDED NOW BECAUSE IT IS NO LONGER AVAILABE FROM www.cleondris.ch 
                else
                    fprintf('Downloading %s\nFrom %s\n',ganymed_java_library_zip,ganymed_java_library_http);
                    urlwrite(ganymed_java_library_http, ganymed_java_library_zip);
                end
                fprintf('Unzipping %s\n',ganymed_java_library_zip);
                unzip(ganymed_java_library_zip,fileparts(ganymed_java_library_zip));
                if (exist(ganymed_java_library_jar_path,'file'))
                    fprintf('Adding Ganymed-ssh2 to the java path by running\njavaaddpath(''%s'')\n',ganymed_java_library_jar_path);
                    javaaddpath(ganymed_java_library_jar_path);
                    fprintf('\nJust added Ganymed-ssh2 to Matlab''s dynamic java Classpath.\n');
                end
            else
                fprintf('A custom ganymed-ssh2 (%s) library should be used and wasn''t, it should be included with the ssh2 matlab files.\n',ganymed_java_library);
                fprintf('Otherwise, change the variable USE_CUSTOM_GANYMED_LIB to zero to automatically download the offical ganymed-ssh2 version');
                fprintf(', which doesn''t support SFTP-get files.');
            end
        end
        [ssh2_struct] = ssh2_config_check_for_java_lib(ssh2_struct);
    end
else
    ssh2_struct = [];
end


function [ssh2_struct] = ssh2_config_check_for_java_lib(ssh2_struct)
jcp = javaclasspath('-all');
gany_found = 0;
for jcp_index = 1:length(jcp)
    if size(strfind(jcp{jcp_index},ssh2_struct.ganymed_java_library_jar_path),1)
        gany_found = gany_found + 1;
    end
end
if (gany_found > 0)
    ssh2_struct.ssh2_java_library_loaded = 1; % lib is loaded
else
    ssh2_struct.ssh2_java_library_loaded = 0; % lib is NOT loaded
end



