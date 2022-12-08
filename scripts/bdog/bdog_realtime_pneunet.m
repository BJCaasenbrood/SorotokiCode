clr;
%% assign free DOF
port = 8888;
usr = 'pi';
pwd = 'softroboticsSA';
ip  = '192.168.0.2';

brd = Bdog(usr,ip,pwd,'autoConnect',false,...
    'Port',port,'NVeab',1);

%% desired pressure
brd.P0  = [2.3,-0.29];
brd.NVeab = 1;
Pd = 30;

%% get data
brd = brd.set('Frequency',25);

%%
pipe = realsense.pipeline();
pointcloud = realsense.pointcloud();
colorizer = realsense.colorizer();
profile = pipe.start();
%dev = profile.get_device();

%% execute control loop
ii = 1;
Node = cell(375,1);
IMG = cell(375,1);

while brd.loop(15)

    fs = pipe.wait_for_frames();
    
    % read data
    %data = brd.tcpRecvData(1);
    t = brd.t;

    p1 = Pd*sigmoid(t)*(0.5+0.5*sin(t));
    %p2 = Pd*sigmoid(t)*sin(t - (2/3)*pi);
    %p3 = Pd*sigmoid(t)*sin(t - (4/3)*pi);
    
    P = [p1,0];
    
    brd.setInput(P);

    depth = fs.get_depth_frame();

    color = colorizer.colorize(depth);

    data = color.get_data();
    %img = permute(reshape(data',[3,color.get_width(),color.get_height()]),[3 2 1]);

    Node{ii} = depth2nodes(pointcloud,depth,[0.35,1,-0.1,0.1,-0.1,.06]);
    img = permute(reshape(data',[3,color.get_width(),color.get_height()]),[3 2 1]);
    img = img(190:350,350:500,:);
    IMG{ii} = img;
    
    fprintf('time = %2.3f\n',brd.t);
    ii = ii + 1;
end

%depth = fs.get_depth_frame();

brd.tcpSendData([0.50, 0.50]); 
%brd.reset();
pause(0.1);

brd.disconnect();
%%
%depth = fs.get_depth_frame();
color = colorizer.colorize(depth);

data = color.get_data();
img = permute(reshape(data',[3,color.get_width(),color.get_height()]),[3 2 1]);
img = img(190:350,350:500,:);

%Node = depth2nodes(pointcloud,depth,[0.35,1,-0.1,0.1,-0.1,.06]);
%Node = depth2nodes(pointcloud,depth);
%Node = depth2nodes(points);

%% plotting
figure(106); clf;
subplot(1,2,1);
t = brd.Log.t;
Y = brd.Log.y(:,1:2);
plot(t,yd(t,Pd)); hold on;
plot(t,brd.Log.y(:,1:2));

%%
for a = 1:375
subplot(1,2,2); 
%scatter3(Node{a}(:,1),Node{a}(:,2),Node{a}(:,3),1,Node{a}(:,1));
imshow(IMG{a});
axis equal;
%view(-90,0);
axis off;
pause(1/25);
drawnow limitrate;
end
% %imshow(img)
%%
filename = ['U_',num2str(Pd),'.mat'];
clearvars -except Node IMG filename t Y
% % 
save(filename);



