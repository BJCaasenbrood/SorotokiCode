%% SSH2 EXAMPLES 
%
% a list of example commands and best practices for the SSH2 functions.
% Please read carefully. Direct any questions to dfreedma@bu.edu.

%% PLEASE SEND ME A DESCRIPTION OF YOUR USE OF SSH2.M. It's quite rewarding
%% and makes the effort and support worthwhile to hear this.



%% IMPORTANT DETAILS
% almost every function returns a ssh2_conn struct. This struct contains
% everything one needs to work with the ssh2 library files. The exception
% to this are the ssh2_simple_command(...), which just returns response from the
% host, the ssh2_command(...)_response, which also retuns the response from the 
% host, and the ssh2_command(...) which returns the ssh2_conn struct and the
% response from the host. ssh2(...) can also return the command response
% when called in a overloaded way (as described near the bottom of this 
% document). 

% ssh_conn - This struct should be reused and saved with all
% of these functions. For more information, please see the non-simple
% examples below.



%% SIMPLE CONNECTION
% Basic Connection: to make a simple, one-time use, connection.
command_output = ssh2_simple_command(HOSTNAME,USERNAME,PASSWORD,'ls -la *ninjas*');
% this command makes a remote ssh connection to the host HOSTNAME, 
% with username USERNAME, and password PASSWORD. 
% Once connected the command 'ls -la *ninjas*' looked
% for ninja files on the remote host. The output of this command will be
% printed in the matlab command window and returned to the variable command_output.


%% SIMPLE CONNECTION WITH OUTPUT TO MATLAB FORCED
% Basic Connection, Output Command Response to MATLAB: If no output arguments
% are specified, the output from the host is displayed in MATLAB. If this is 
% desired, just provide an additional input set to 1.
command_output = ssh2_simple_command(HOSTNAME,USERNAME,PASSWORD,'ls -la *ninjas*',1);


%% ADVANCED CONNECTIONS
% The simple commands will connect, authenticate, issue the command, and
% then close the connection. For those who need to issue multiple commands
% or need more advanced options, it's recommend to setup a ssh2_conn
% configuration and use the non-simple ssh2 functions. 
% For example, to perform the same thing as before, but have the ssh2_conn
% configuration available for multiple connections. It is best practice to
% return the ssh2_conn to the same variable. This way the connection
% information is always current and correct. 
% Advanced Connection: 
ssh2_conn = ssh2_config(HOSTNAME,USERNAME,PASSWORD);
ssh2_conn = ssh2_command(ssh2_conn, 'ls -la *ninjas*');
% or to suppress the output from being displayed I can run
% Reusing Advanced Connection: 
ssh2_conn = ssh2_command(ssh2_conn, 'ls -la *ninjas*',1);
% the output is available in cell array ssh2_conn.command_result.
% if I want to run multiple commands, I can just string them together with 
% the command seperator, the semicolon
% RETRIEVING THE COMMAND RESPONSE:
% the command response can be accessed easily through the
% ssh2_command_response() method. 
command_response = ssh2_command_response(ssh2_conn);
% to access the first response, one can use:
command_response{1}
% to get the number of lines returned
size(command_response,1);


% Multiple Commands, Advanced Connection: 
ssh2_conn = ssh2_command(ssh2_conn, 'ls *dogs*; ls *cats*',1);
% when I'm done, I need to make sure to close an advanced connection
% Close Advanced Connection: 
ssh2_conn = ssh2_close(ssh2_conn);
% NOTE: Most values for the command, scp, and sftp are cleared after 
% each use. This is OK as long as one doesn't then manually change/or not
% change the values of the ssh2_conn struct and call ssh2 directly. 


%% PUBLIC/PRIVATE KEY AUTHENTICATION
% If one would prefer to use public/private keys for authentication, just
% use the advnaced connection and generate your ssh2_conn config using the
% ssh2_config_publickey function instead. Everything else will be
% identical.
%%% !!! You have to convert the key into openssh format (*.pem file)

ssh2_conn = ssh2_config_publickey(HOSTNAME,USERNAME,PRIVATE_KEY_FILE,PRIVATE_KEY_PASSWORD);
[ssh2_conn,command_output] = ssh2_command(ssh2_conn, 'ls -la *birds*'); %get output in second variable
ssh2_conn = ssh2_command(ssh2_conn, 'ls -la *bess*',1); %suppress outpiut
ssh2_conn = scp_get(ssh2_conn, 'readme.txt'); %scp file from remote host
ssh2_conn = ssh2_close(ssh2_conn); %close connection when done



%% SCP/SFTP
% speaking of file transers, one has a few options. The simple file
% transers functions, scp_simple_get, scp_simple_put, sftp_simple_get,
% sftp_simple_put work the same way as the ssh2_simple_command. One
% supplies the full details of the connection and a one-time connection is
% made.

% SCP is the recommended interface as it is ususally faster and
% because SFTP GET is not supported by these functions, unless the 
% variable USE_CUSTOM_GANYMED_LIB = 1 in the beginning of ssh2_setup.m.
% if USE_CUSTOM_GANYMED_LIB = 1, an included custom ganymed ssh2 library
% (ganymed-ssh2-m1) will be used which does support SFTP_GET.


%% SCP SIMPLE GET
% download a file, "filetodownload.txt" from remote host HOSTNAME
ssh2_conn = scp_simple_get(HOSTNAME,USERNAME,PASSWORD,'filetodownload.txt');
% if one would like to download a few files sequentilly, use a cell array
ssh2_conn = scp_simple_get(HOSTNAME,USERNAME,PASSWORD,{'file1','file2','file3'});
% if one needs to specify a directory to download besides the working
% directory add a fifth input to the function
ssh2_conn = scp_simple_get(HOSTNAME,USERNAME,PASSWORD,'filetodownload.txt','/path/to/download/to');
% or if the remote file is in a directory outside of home
ssh2_conn = scp_simple_get(HOSTNAME,USERNAME,PASSWORD,{'file1','file2','file3'},'/path/to/download/to','remotehostpath/folder');


%% SCP SIMPLE PUT
% upload a file, "filetoupload.txt" to from remote host HOSTNAME
ssh2_conn = scp_simple_put(HOSTNAME,USERNAME,PASSWORD,'filetoupload.txt');
% if one would like to upload a few files sequentilly, use a cell array
ssh2_conn = scp_simple_put(HOSTNAME,USERNAME,PASSWORD,{'file1','file2','file3'});
% if one needs to specify a directory to upload to besides the default directory 
ssh2_conn = scp_simple_put(HOSTNAME,USERNAME,PASSWORD,'filetoupload.txt','/path/to/upload/to');
% or if the local file is in a directory outside of the working directory
ssh2_conn = scp_simple_put(HOSTNAME,USERNAME,PASSWORD,'filetoupload.txt','/path/to/upload/to','C:\localfolder\');
% and finally, if one would like to rename the file when uploading, supply
% an additional input value
ssh2_conn = scp_simple_put(HOSTNAME,USERNAME,PASSWORD,'filetoupload.txt','/path/to/upload/to','C:\localfolder\','uploadedfile_renamed.txt');
% if using a cell array to upload files and you wish to rename them, it
% should be pretty obvious you need to use a cell array for the renamed
% files. If you don't, I have no idea what will happen... And why?
ssh2_conn = scp_simple_put(HOSTNAME,USERNAME,PASSWORD,{'file1','file2','file3'},'/path/to/upload/to','C:\localfolder\',{'apple','banana','orange'});


%% SFTP
%  anywhere scp is used, sftp can be replaced. The input syntax is the same.
% NOTE: SCP is the recommended interface as it is ususally faster and
% because SFTP GET is not supported by these functions, unless the 
% variable USE_CUSTOM_GANYMED_LIB = 1 in set in the beginning of ssh2_setup.m.
% if USE_CUSTOM_GANYMED_LIB = 1, an included custom ganymed ssh2 library
% (ganymed-ssh2-m1) will be used which does support SFTP_GET.

%% ADVANCED SCP GET/PUT
% advanced connections that authenticate once and can be reused are
% recommneded for more than one connection. These are very similar to the
% advanced command connections. The last three (get) and four (put)
% variables in these functions are idenctical the simple cases explained
% above.
ssh2_conn = ssh2_config(HOSTNAME,USERNAME,PASSWORD);
ssh2_conn = ssh2_command(ssh2_conn, 'ls -la *ninjas*');
ssh2_conn = scp_get(ssh2_conn, {'file1','file2','file3'},'/path/to/download/to','remotehostpath/folder');
ssh2_conn = scp_put(ssh2_conn, 'filetoupload.txt','/path/to/upload/to','C:\localfolder\','uploadedfile_renamed.txt');

%% ADVANCED CLOSING A CONNECTION
% advanced close below. This illustrates how an advanced connection be
% manipulated and how ssh2_close calls ssh2 and then closes the connection.
% This is the most optimized way to close the conneciton because it doesn't
% require checking the input variables of ssh2 again (when called in
% ssh2_command).
ssh2_conn.command = 'ls -la filetoupload.txt'; % checking that new file is there
ssh2_conn = ssh2_close(ssh2_conn); %will call ssh2.m and run command and then close connection
% commands are issued after scp and sftp. This would be a perfect way to
% upload and check that the file was uploaded sucessfully, or perform an
% operation on the file (chmod +x)

%% ADVANCED, UNDOCUMENTED FEATURES
% there are a few options that are not directly marketed

% PORT, non-standard ssh port
% Add a fourth variable ssh_config, or a fifth to ssh2_config_publickey
ssh2_conn = ssh2_config(HOSTNAME,USERNAME,PASSWORD,PORT);
ssh2_conn = sh2_config_publickey(HOSTNAME,USERNAME,PRIVATE_KEY_FILE,PRIVATE_KEY_PASSWORD,PORT);
% also after a ssh2_conn config has been created, one can modify it
% directly
ssh2_conn.port = 2200; %non-standard ssh port

% COMMAND RESPONSE
% the command response from the host can be found in the cell array
ssh2_conn.command_result

% AUTO-RECONNECT
% there is an auto-reconnect feature. By settting autoreconnect = 1, ssh2
% will check if it's connected everytime it is called. This will require
% a session to be created and closed everytime, before any other function
% is done. It will then attenpt to reconnect if a connection is somehow
% broken. Can be useful, but adds overhead. Disabled by default.
ssh2_conn.autoreconnect = 1; % to enable auto-reconnect option

%% SCP MODE
% when using SCP to copy files, it is possible to establish the permissions
% of newly created files. This can be accomplished by setting
% remote_file_mode to an integer mode value, that is 0644 for rw-r--r--
% permissions. Default is 0600. This is not used with SFTP.
ssh2_conn.remote_file_mode = 0644; % change to a more readable permission

%% IGNORE HOST RESPONSE
% for the fastest connections, one can send the commands and not readback
% the host response. 
ssh2_conn.command_ignore_response = 1;

% it should be stated that ssh2(...) (or ssh2_main(...), but this should be called 
% through ssh2.m) is the bulk of the ssh2 engine. ssh2.m is NOT intended 
% to be called directly, although this is possible with a correctly 
% configured ssh2_conn struct. For convenience, ssh2(...) can be used
% as with the same inputs as ssh2_simple_command and will run this command
% instead, e.g.
ssh2_conn = ssh2(HOSTNAME,USERNAME,PASSWORD,'ls -la *.m');
% will call ssh2_conn = ssh2_simple_command(HOSTNAME,USERNAME,PASSWORD,'ls -la *.m');
%
% one can also use another ssh2 shorthand. If one has a configured
% ssh2_conn, provide that and a command to run, as the second argument, and
% one can run ssh2_command directly from ssh2, e.g.
ssh2_conn = ssh2(ssh2_conn,'ls -la *.m');
% will call ssh2_conn = ssh2_command(ssh2_conn,'ls -la *.m'); This
% invocation will never output the command response to the screen.


% ssh2_main(...) is NEVER meant to be called directly. Always call this function
% through ssh2(...)