clr;
%% assign free DOF
port = 8888;
usr = 'pi';
pwd = 'softroboticsSA';
ip  = '192.168.0.2';

brd = Bdog(usr,ip,pwd,'autoConnect',false,'Port',port);


%% desired pressure
brd.NVeab = 2;
Pd = 10;

%% get data
brd = brd.set('Frequency',50);

%% execute control loop

while brd.loop(10)
    
    % read data
    data = brd.tcpRecvData(1);
    
    Pmeasure = padc(data(1),0.4);

    % control law: 
    u = ControlLaw(brd.t,Pd,Pmeasure);
    
    % send data 
    brd.tcpSendData([u, 0.5, 0.5, 0.5]);
    
    fprintf('time = %2.3f\n',brd.t);
end

brd.tcpSendData([0.5, 0.5, 0.45, 0.5]); pause(0.1);

brd.disconnect();

%% plotting
t = brd.Log.t;
plot(t,yd(t,Pd)); hold on;
plot(t,padc(brd.Log.y(:,1),0.4));

%% additional functions
function u = ControlLaw(t,Pd,Pm)
    
    % P-gain
    Kp = 0.25;
    
    % pressure-to-voltage
    V = pdac(yd(t,Pd));
    
    % control action
    u = V - Kp*(pdac(Pm)-V);
end

function z = yd(t,Pd)
    z = Pd*sin(t);%lerp(0,Pd,smoothstep(0.75*t));
end


