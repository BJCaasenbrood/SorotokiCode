clr;
%% assign free DOF
port = 8888;
usr = 'pi';
pwd = 'softroboticsSA';
ip  = '192.168.0.2';

brd = Bdog(usr,ip,pwd,'autoConnect',false,...
    'Port',port,'NVeab',2);

%% desired pressure
brd.P0  = [1,1,0,0];
brd.NVeab = 2;
Pd = 10;

%% get data
brd = brd.set('Frequency',50);

%% execute control loop

while brd.loop(15)
    
    % read data
    %data = brd.tcpRecvData(1);
    t = brd.t;
    p1 = Pd*sigmoid(t)*sin(t - (0/3)*pi);
    p2 = Pd*sigmoid(t)*sin(t - (2/3)*pi);
    p3 = Pd*sigmoid(t)*sin(t - (4/3)*pi);
    
    P = [p1,p2,p3,0];
    %brd.setInput([p1,p2,p3,0]);
    brd.setInput(P+[3,3,0.5,0.5]);
    %P = [p1,p2,p3,0].';
    
    %brd.setInput(P);
    
    %Pmeasure = brd.getOutput();
    
%     Pd = [0,0,0,0];
%     %U = p
% %     Pmeasure = padc(data(1),0.4);
% % 
% %     % control law: 
     %u = ControlLaw(brd.t,Pd,Pmeasure);
%     
%     % send data 
%     brd.tcpSendData([1,1,1,1]*0.5);
    
    fprintf('time = %2.3f\n',brd.t);
end

brd.tcpSendData([0.508, 0.508, 0.5, 0.5]); 
%brd.reset();
pause(0.1);

brd.disconnect();

%% plotting
t = brd.Log.t;
clf;
%plot(t,yd(t,Pd)); hold on;
plot(t,padc(brd.Log.y(:,1:3),0.4));

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
    z = Pd;%*sin(t);
end


