function SSH_ = connectBoard(Bdog)

ip = Bdog.ip;   %'169.254.168.25';
usr = Bdog.usr; %'pi';
pwd = Bdog.pwd; %'sorotoki'
fprintf('~finding host computer... \n'); 

SSH_ = ssh2_config(ip,usr,pwd);
SSH_ = ssh2(SSH_);

cout('green',['* Connected to host computer: @',ip,'\n']);
end