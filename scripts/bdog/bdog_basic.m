clr;
%% assign free DOF
usr = 'pi';
pwd = 'softroboticsSA';
ip  = '192.168.0.2';

brd = Bdog(usr,ip,pwd,'autoConnect',true,'Port',8889);

%% get data
brd = brd.set('Frequency',150);

%%
%disp('RUN: "python3 main.py 3" on the RPI');
%brd = brd.command('cd Sorotoki && nohup python3 main.py');
pause(0.5);

%% execute control loop
while brd.loop(10)
    
    % read data
    data = brd.tcpRecvData(1);
    
    % desired pressure in kPa
    pd = 20*(0.5*(sin(brd.t) + 1));
    
    % pressure to Digital-Analog voltage
    V = pdac(pd);
    
    % send data 
    brd.tcpSendData([0.5, V]);
    
end

% reset VEAB to 0.5V
brd.tcpSendData([0.5, 0.5]);
brd.disconnect();

%% plotting
t = brd.Log.t;
y = brd.Log.y;
%plot(t,y); hold on;
% 
% ADC to pressure measurements
Pm = padc(y);

plot(t,Pm);

