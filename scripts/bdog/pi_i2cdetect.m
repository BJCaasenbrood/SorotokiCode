clr;
%% assign free DOF
usr = 'pi';
pwd = 'softroboticsSA';
ip  = '192.168.0.2';

%% 
brd = Bdog(usr,ip,pwd,'autoConnect',true);

%% loop of I2C ports
cmd = ['cd Documents/Sorotoki && python3 i2cfind.py']; 

brd = brd.command(cmd);   

brd.disconnect();