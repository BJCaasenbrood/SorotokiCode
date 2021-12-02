clr;
%% assign free DOF
usr = 'pi';
pwd = 'softroboticsSA';
ip  = '192.168.0.2';

brd = Bdog(usr,ip,pwd,'autoConnect',true);

%% get data
brd = brd.set('Frequency',400);
% brd.transfer('baseSoftRobot.py');
% brd.transfer('SoftRobot.py');
% brd.transfer('runme.py');

%% execute control loop
Kp = 0.65;

while brd.loop(5)
    
    % read data
    data = brd.tcpRecvData(1);

    % control law: 
    %u = yd(brd.t) - Kp*(data/10-yd(brd.t));
    u = yd(brd.t);
    
    % send data 
    brd.tcpSendData(u);
end

brd.disconnect();

%% plotting
% t = brd.Log.Time;
% plot(t,yd(t)); hold on;
% plot(brd.Log.Time,brd.Log.Data/10);

%%dy additional functions
function y = yd(t)
    y = 0.4 + 0.1*sin(t);
end