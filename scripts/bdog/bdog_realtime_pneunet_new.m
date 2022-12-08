clr;
%% assign free DOF
port = 8888;
usr = 'pi';
pwd = 'softroboticsSA';
ip  = '192.168.0.2';

brd = Bdog(usr,ip,pwd,'Port',port,'NVeab',1);

%% camera
cam = Vision('realsense');
cam = cam.set('HoodXY',21,'HoodR',51,...
    'Frame',[190,553,137,409],...
    'ColorFilter',[10,40],...
    'UpdateGain',6/25);

cam = cam.setWorld([0,260]-107.55,[0,195]-55.35);

%% desired pressure
brd.P0    = [2.3,-0.29];
brd.NVeab = 1;
Pd = 30;

%% get data
brd = brd.set('Frequency',25);

%% execute control loop
h = [];
Y = [];

while brd.loop(20)
    
    % read data
    t = brd.t;

    p1 = Pd*sigmoid(t)*(0.5+0.5*sin(t));
    
    P = [p1,0];
    
    brd.setInput(P);

    cam = cam.getframe();

    if isempty(h)
        h = imshow(cam.lastFrame);
    else
        h.CData = cam.lastFrame;
    end

    cam = cam.detect(9);
    Y = [Y;cam.output(:).'];
    
    fprintf('time = %2.3f\n',brd.t);
    
end

brd.tcpSendData([0.50, 0.50]); 

pause(0.1);

brd.disconnect();
%%
%Node = depth2nodes(pointcloud,depth,[0.35,1,-0.1,0.1,-0.1,.06]);
%Node = depth2nodes(pointcloud,depth);
%Node = depth2nodes(points);

figure(103);
plot(Y(:,1),'LineWidth',1.5); hold on;
plot(Y(:,2),'LineWidth',1.5); hold on;
sorocolor;





