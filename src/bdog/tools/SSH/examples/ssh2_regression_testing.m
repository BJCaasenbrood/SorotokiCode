function ssh2_regression_testing()
%% This file is not 100% complete
HOSTNAME = '192.168.0.2';
USERNAME = 'username';
PASSWORD = 'password';
cmd = 'ls -la ssh2.zip';
local_file = 'ganymed-ssh2-build250.zip';
local_file_folder = 'PhotomaskBasics-General.pdf';
local_folder = '/Users/david';
remote_file = 'ByteArrayBuffer.jar';
remote_file_folder = 'ssh2.zip';
remote_folder = '/dump/dev/matlab';
remote_file_new_name = 'gany.zip';
PEMFILE = '';
PEMFILEPASSWORD = '';

passed = 0;

%% SSH2 SETUP
conn = ssh2_setup();
conn = ssh2_setup(conn);
conn = ssh2_setup(now());

%% BASIC CONNECTION
% invalid input
command_output = ssh2_simple_command(HOSTNAME,USERNAME,PASSWORD);
% correct connection
command_output = ssh2_simple_command(HOSTNAME,USERNAME,PASSWORD,cmd);
% no output to be displayed
ssh2_simple_command(HOSTNAME,USERNAME,PASSWORD,cmd);
command_output = ssh2_simple_command(HOSTNAME,USERNAME,PASSWORD,cmd,1);
% no output to be displayed
ssh2_simple_command(HOSTNAME,USERNAME,PASSWORD,cmd,1);

%%SSH2 CONFIG
ssh2conn = ssh2_config_publickey(HOSTNAME,USERNAME);
ssh2conn = ssh2_config_publickey(HOSTNAME,USERNAME,PEMFILE,PEMFILEPASSWORD);
ssh2conn = ssh2_close();
ssh2conn = ssh2_close(ssh2conn);

%%SSH2 CONFIG
ssh2conn = ssh2_config(HOSTNAME,USERNAME);
ssh2conn = ssh2_config(HOSTNAME,USERNAME,PASSWORD);
ssh2conn = ssh2_close(ssh2conn);

ssh2_conn = ssh2_config(HOSTNAME, USERNAME, PASSWORD);

extrastr = '';
ssh2conn = sftp_simple_put(HOSTNAME, USERNAME, PASSWORD);
[filethere] = check_remote_file(ssh2_conn, local_file,1); fprintf('File Threre: %d\n',filethere);
ssh2conn = sftp_simple_put(HOSTNAME, USERNAME, PASSWORD, local_file);
[filethere] = check_remote_file(ssh2_conn, local_file,1); fprintf('File Threre: %d\n',filethere);
ssh2conn = sftp_simple_put(HOSTNAME, USERNAME, PASSWORD, local_file, [remote_folder extrastr]);
[filethere] = check_remote_file(ssh2_conn, [remote_folder extrastr local_file],1); fprintf('File Threre: %d\n',filethere);
ssh2conn = sftp_simple_put(HOSTNAME, USERNAME, PASSWORD, local_file_folder, [remote_folder extrastr], [local_folder extrastr]);
[filethere] = check_remote_file(ssh2_conn, [remote_folder extrastr local_file_folder],1); fprintf('File Threre: %d\n',filethere);
ssh2conn = sftp_simple_put(HOSTNAME, USERNAME, PASSWORD, local_file_folder, [remote_folder extrastr], [local_folder extrastr], remote_file_new_name);
[filethere] = check_remote_file(ssh2_conn, [remote_folder extrastr remote_file_new_name],1); fprintf('File Threre: %d\n',filethere);

extrastr = '/';
ssh2conn = sftp_simple_put(HOSTNAME, USERNAME, PASSWORD);
[filethere] = check_remote_file(ssh2_conn, local_file,1); fprintf('File Threre: %d\n',filethere);
ssh2conn = sftp_simple_put(HOSTNAME, USERNAME, PASSWORD, local_file);
[filethere] = check_remote_file(ssh2_conn, local_file,1); fprintf('File Threre: %d\n',filethere);
ssh2conn = sftp_simple_put(HOSTNAME, USERNAME, PASSWORD, local_file, [remote_folder extrastr]);
[filethere] = check_remote_file(ssh2_conn, [remote_folder extrastr local_file],1); fprintf('File Threre: %d\n',filethere);
ssh2conn = sftp_simple_put(HOSTNAME, USERNAME, PASSWORD, local_file_folder, [remote_folder extrastr], [local_folder extrastr]);
[filethere] = check_remote_file(ssh2_conn, [remote_folder extrastr local_file_folder],1); fprintf('File Threre: %d\n',filethere);
ssh2conn = sftp_simple_put(HOSTNAME, USERNAME, PASSWORD, local_file_folder, [remote_folder extrastr], [local_folder extrastr], remote_file_new_name);
[filethere] = check_remote_file(ssh2_conn, [remote_folder extrastr remote_file_new_name],1); fprintf('File Threre: %d\n',filethere);


extrastr = '';
ssh2conn = scp_simple_put(HOSTNAME, USERNAME, PASSWORD);
[filethere] = check_remote_file(ssh2_conn, local_file,1); fprintf('File Threre: %d\n',filethere);
ssh2conn = scp_simple_put(HOSTNAME, USERNAME, PASSWORD, local_file);
[filethere] = check_remote_file(ssh2_conn, local_file,1); fprintf('File Threre: %d\n',filethere);
ssh2conn = scp_simple_put(HOSTNAME, USERNAME, PASSWORD, local_file, [remote_folder extrastr]);
[filethere] = check_remote_file(ssh2_conn, [remote_folder extrastr local_file],1); fprintf('File Threre: %d\n',filethere);
ssh2conn = scp_simple_put(HOSTNAME, USERNAME, PASSWORD, local_file_folder, [remote_folder extrastr], [local_folder extrastr]);
[filethere] = check_remote_file(ssh2_conn, [remote_folder extrastr local_file_folder],1); fprintf('File Threre: %d\n',filethere);
ssh2conn = scp_simple_put(HOSTNAME, USERNAME, PASSWORD, local_file_folder, [remote_folder extrastr], [local_folder extrastr], remote_file_new_name);
[filethere] = check_remote_file(ssh2_conn, [remote_folder extrastr remote_file_new_name],1); fprintf('File Threre: %d\n',filethere);

extrastr = '/';
ssh2conn = scp_simple_put(HOSTNAME, USERNAME, PASSWORD);
[filethere] = check_remote_file(ssh2_conn, local_file,1); fprintf('File Threre: %d\n',filethere);
ssh2conn = scp_simple_put(HOSTNAME, USERNAME, PASSWORD, local_file);
[filethere] = check_remote_file(ssh2_conn, local_file,1); fprintf('File Threre: %d\n',filethere);
ssh2conn = scp_simple_put(HOSTNAME, USERNAME, PASSWORD, local_file, [remote_folder extrastr]);
[filethere] = check_remote_file(ssh2_conn, [remote_folder extrastr local_file],1); fprintf('File Threre: %d\n',filethere);
ssh2conn = scp_simple_put(HOSTNAME, USERNAME, PASSWORD, local_file_folder, [remote_folder extrastr], [local_folder extrastr]);
[filethere] = check_remote_file(ssh2_conn, [remote_folder extrastr local_file_folder],1); fprintf('File Threre: %d\n',filethere);
ssh2conn = scp_simple_put(HOSTNAME, USERNAME, PASSWORD, local_file_folder, [remote_folder extrastr], [local_folder extrastr], remote_file_new_name);
[filethere] = check_remote_file(ssh2_conn, [remote_folder extrastr remote_file_new_name],1); fprintf('File Threre: %d\n',filethere);

extrastr = '';
ssh2conn = scp_put(ssh2conn);
[filethere] = check_remote_file(ssh2_conn, local_file,1); fprintf('File Threre: %d\n',filethere);
ssh2conn = scp_put(ssh2conn, local_file);
[filethere] = check_remote_file(ssh2_conn, local_file,1); fprintf('File Threre: %d\n',filethere);
ssh2conn = scp_put(ssh2conn, local_file, [remote_folder extrastr]);
[filethere] = check_remote_file(ssh2_conn, [remote_folder extrastr local_file],1); fprintf('File Threre: %d\n',filethere);
ssh2conn = scp_put(ssh2conn, local_file_folder, [remote_folder extrastr], [local_folder extrastr]);
[filethere] = check_remote_file(ssh2_conn, [remote_folder extrastr local_file_folder],1); fprintf('File Threre: %d\n',filethere);
ssh2conn = scp_put(ssh2conn, local_file_folder, [remote_folder extrastr], [local_folder extrastr], remote_file_new_name);
[filethere] = check_remote_file(ssh2_conn, [remote_folder extrastr remote_file_new_name],1); fprintf('File Threre: %d\n',filethere);

extrastr = '/';
ssh2conn = scp_put(ssh2conn);
[filethere] = check_remote_file(ssh2_conn, local_file,1); fprintf('File Threre: %d\n',filethere);
ssh2conn = scp_put(ssh2conn, local_file);
[filethere] = check_remote_file(ssh2_conn, local_file,1); fprintf('File Threre: %d\n',filethere);
ssh2conn = scp_put(ssh2conn, local_file, [remote_folder extrastr]);
[filethere] = check_remote_file(ssh2_conn, [remote_folder extrastr local_file],1); fprintf('File Threre: %d\n',filethere);
ssh2conn = scp_put(ssh2conn, local_file_folder, [remote_folder extrastr], [local_folder extrastr]);
[filethere] = check_remote_file(ssh2_conn, [remote_folder extrastr local_file_folder],1); fprintf('File Threre: %d\n',filethere);
ssh2conn = scp_put(ssh2conn, local_file_folder, [remote_folder extrastr], [local_folder extrastr], remote_file_new_name);
[filethere] = check_remote_file(ssh2_conn, [remote_folder extrastr remote_file_new_name],1); fprintf('File Threre: %d\n',filethere);

extrastr = '';
ssh2conn = sftp_put(ssh2conn);
[filethere] = check_remote_file(ssh2_conn, local_file,1); fprintf('File Threre: %d\n',filethere);
ssh2conn = sftp_put(ssh2conn, local_file);
[filethere] = check_remote_file(ssh2_conn, local_file,1); fprintf('File Threre: %d\n',filethere);
ssh2conn = sftp_put(ssh2conn, local_file, [remote_folder extrastr]);
[filethere] = check_remote_file(ssh2_conn, [remote_folder extrastr local_file],1); fprintf('File Threre: %d\n',filethere);
ssh2conn = sftp_put(ssh2conn, local_file_folder, [remote_folder extrastr], [local_folder extrastr]);
[filethere] = check_remote_file(ssh2_conn, [remote_folder extrastr local_file_folder],1); fprintf('File Threre: %d\n',filethere);
ssh2conn = sftp_put(ssh2conn, local_file_folder, [remote_folder extrastr], [local_folder extrastr], remote_file_new_name);
[filethere] = check_remote_file(ssh2_conn, [remote_folder extrastr remote_file_new_name],1); fprintf('File Threre: %d\n',filethere);

extrastr = '/';
ssh2conn = sftp_put(ssh2conn);
[filethere] = check_remote_file(ssh2_conn, local_file,1); fprintf('File Threre: %d\n',filethere);
ssh2conn = sftp_put(ssh2conn, local_file);
[filethere] = check_remote_file(ssh2_conn, local_file,1); fprintf('File Threre: %d\n',filethere);
ssh2conn = sftp_put(ssh2conn, local_file, [remote_folder extrastr]);
[filethere] = check_remote_file(ssh2_conn, [remote_folder extrastr local_file],1); fprintf('File Threre: %d\n',filethere);
ssh2conn = sftp_put(ssh2conn, local_file_folder, [remote_folder extrastr], [local_folder extrastr]);
[filethere] = check_remote_file(ssh2_conn, [remote_folder extrastr local_file_folder],1); fprintf('File Threre: %d\n',filethere);
ssh2conn = sftp_put(ssh2conn, local_file_folder, [remote_folder extrastr], [local_folder extrastr], remote_file_new_name);
[filethere] = check_remote_file(ssh2_conn, [remote_folder extrastr remote_file_new_name],1); fprintf('File Threre: %d\n',filethere);



%%GET 
extrastr = '';
ssh2conn = sftp_simple_get(HOSTNAME, USERNAME, PASSWORD);
[filethere] = check_local_file(local_file,1); fprintf('File Threre: %d\n',filethere);
ssh2conn = sftp_simple_get(HOSTNAME, USERNAME, PASSWORD, remote_file);
[filethere] = check_local_file(local_file,1); fprintf('File Threre: %d\n',filethere);
ssh2conn = sftp_simple_get(HOSTNAME, USERNAME, PASSWORD, remote_file, [local_folder extrastr]);
[filethere] = check_local_file([local_folder extrastr local_file],1); fprintf('File Threre: %d\n',filethere);
ssh2conn = sftp_simple_get(HOSTNAME, USERNAME, PASSWORD, remote_file, [local_folder extrastr], [remote_folder extrastr]);
[filethere] = check_local_file([local_folder extrastr local_file_folder],1); fprintf('File Threre: %d\n',filethere);


extrastr = '/';




return;

command_output = ssh2_simple_command(HOSTNAME,USERNAME,PASSWORD,'ls -la *ninjas*',1);

ssh2_conn = ssh2_config(HOSTNAME,USERNAME,PASSWORD);
ssh2_conn = ssh2_command(ssh2_conn, 'ls -la *ninjas*');
ssh2_conn = ssh2_command(ssh2_conn, 'ls -la *ninjas*',1);
command_response = ssh2_command_response(ssh2_conn);
ssh2_conn = ssh2_command(ssh2_conn, 'ls *dogs*; ls *cats*',1);
ssh2_conn = ssh2_close(ssh2_conn);
ssh2_conn = ssh2_config_publickey(HOSTNAME,USERNAME,PRIVATE_KEY_FILE,PRIVATE_KEY_PASSWORD);
[ssh2_conn,command_output] = ssh2_command(ssh2_conn, 'ls -la *birds*'); %get output in second variable
ssh2_conn = ssh2_command(ssh2_conn, 'ls -la *bess*',1); %suppress outpiut
ssh2_conn = scp_get(ssh2_conn, 'readme.txt'); %scp file from remote host
ssh2_conn = ssh2_close(ssh2_conn); %close connection when done
ssh2_conn = scp_simple_get(HOSTNAME,USERNAME,PASSWORD,'filetodownload.txt');
ssh2_conn = scp_simple_get(HOSTNAME,USERNAME,PASSWORD,{'file1','file2','file3'});
ssh2_conn = scp_simple_get(HOSTNAME,USERNAME,PASSWORD,'filetodownload.txt','/path/to/download/to');
ssh2_conn = scp_simple_get(HOSTNAME,USERNAME,PASSWORD,{'file1','file2','file3'},'/path/to/download/to','remotehostpath/folder');

ssh2_conn = scp_simple_put(HOSTNAME,USERNAME,PASSWORD,'filetoupload.txt');
ssh2_conn = scp_simple_put(HOSTNAME,USERNAME,PASSWORD,{'file1','file2','file3'});
ssh2_conn = scp_simple_put(HOSTNAME,USERNAME,PASSWORD,'filetoupload.txt','/path/to/upload/to');
ssh2_conn = scp_simple_put(HOSTNAME,USERNAME,PASSWORD,'filetoupload.txt','/path/to/upload/to','C:\localfolder\');
ssh2_conn = scp_simple_put(HOSTNAME,USERNAME,PASSWORD,'filetoupload.txt','/path/to/upload/to','C:\localfolder\','uploadedfile_renamed.txt');
ssh2_conn = scp_simple_put(HOSTNAME,USERNAME,PASSWORD,{'file1','file2','file3'},'/path/to/upload/to','C:\localfolder\',{'apple','banana','orange'});

ssh2_conn = ssh2_config(HOSTNAME,USERNAME,PASSWORD);
ssh2_conn = ssh2_command(ssh2_conn, 'ls -la *ninjas*');
ssh2_conn = scp_get(ssh2_conn, {'file1','file2','file3'},'/path/to/download/to','remotehostpath/folder');
ssh2_conn = scp_put(ssh2_conn, 'filetoupload.txt','/path/to/upload/to','C:\localfolder\','uploadedfile_renamed.txt');
ssh2_conn.command = 'ls -la filetoupload.txt'; % checking that new file is there
ssh2_conn = ssh2_close(ssh2_conn); %will call ssh2.m and run command and then close connection
% commands are issued after scp and sftp. This would be a perfect way to
ssh2_conn = ssh2_config(HOSTNAME,USERNAME,PASSWORD,PORT);
ssh2_conn = sh2_config_publickey(HOSTNAME,USERNAME,PRIVATE_KEY_FILE,PRIVATE_KEY_PASSWORD,PORT);
ssh2_conn.port = 2200; %non-standard ssh port
ssh2_conn.command_result
ssh2_conn.autoreconnect = 1; % to enable auto-reconnect option
ssh2_conn.remote_file_mode = 0644; % change to a more readable permission
ssh2_conn.command_ignore_response = 1;
ssh2_conn = ssh2(HOSTNAME,USERNAME,PASSWORD,'ls -la *.m');
ssh2_conn = ssh2(ssh2_conn,'ls -la *.m');

function [filethere] = check_remote_file(ssh2_conn, fn,del)
[ssh2_conn, response] = ssh2_command(ssh2_conn,['[ -f ' fn ' ] && echo 1 || echo 0']);
if response{1} == '1'
    filethere = 1;
else
    filethere = 0;
end
if del == 1
   [ssh2_conn, response] = ssh2_command(ssh2_conn,['rm ' fn ';']);
end
 
function [filethere] = check_local_file(fn,del)
if (exist(fn,'file'))
    filethere = 1;
else
    filethere = 0;
end
if del == 1
   delete(fn);
end
 
