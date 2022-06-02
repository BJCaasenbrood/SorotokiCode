clr;
%% assign free DOF
usr = 'pi';
pwd = 'softroboticsSA';
ip  = '192.168.0.2';

brd = Bdog(usr,ip,pwd,'autoConnect',true);

%% get data
brd = brd.set('Frequency',50);

%% execute control loop
while brd.loop(10)
    
    % read data
    data = brd.tcpRecvData(1);
    
    % desired pressure in kPa
    pd = 20*(0.5*(sin(brd.t) + 1));
    
    % pressure to Digital-Analog voltage
    V = pdac(pd);
    
    % send data 
    brd.tcpSendData([V, 0.5]);
    
    % print time
    fprintf('time = %2.3f\n',brd.t);
end

% reset VEAB to 0.5V
brd.tcpSendData([0.5, 0.5]);
brd.disconnect();

%% plotting
t = brd.Log.Time;
plot(t,yd(t,Pd)); hold on;

% ADC to pressure measurements
Pm = padc(brd.Log.Data(:,1));

plot(t,Pm);

