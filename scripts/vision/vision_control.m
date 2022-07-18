clr;
%% build Vision class
figure(101);
cam = Vision(1);

%% seup control-board
ip  = '192.168.0.2';
brd = Bdog('pi',ip,'softroboticsSA',...
    'autoConnect',true);

%% AR detection
cam = cam.detectAR();

%% get data
brd = brd.set('Frequency',10);
shp = SetupModel();

%% detect AR markers
h1   = cam.show(); 
out  = cam.detect(); [x, y] = circlepoints(out(3));
Ctr  = CenterMarker(cam);

%% initial curve
Q0 = [1e-7;1e-7];
Pix = FK2Pix(Q0,Ctr,cam,shp);

%% desired tracjectory
tt   = linspace(0,2*pi,200).';
PixD = Point2Pix([50*sin(0.5*tt),69 + 14*cos(tt)],Ctr,cam);
% 
% PixDc =@(t) Point2Pix([50*sin(0.5*t),69 + 14*cos(t)],Ctr,cam);
% plot(PixD(:,1),PixD(:,2),'k--');

%% input matrix
a = 1;
b = 1;
G = [-a,a;b,b];

%% execute control loop
data = [];  
Pd   = [];
Yd   = [];

hold on;
h2  = plot(x + out(1),y + out(2),'g-','LineW',3);
h22 = plot(out(1),out(2),'go','LineW',3);
h3  = plot(Pix(:,1),Pix(:,2),'o-','LineW',2,'Color',col(1));
drawnow;
h4 = [];

while brd.loop(10)
    
    % read data
    brd.tcpRecvData(1);

    % update cam
    I = snapshot(cam.Cam);
    set(h1,'CData',I);
    
    % control law:     
    out = cam.detect();
    set(h2,'XData',x+out(1),'YData',y+out(2));
    set(h22,'XData',out(1),'YData',out(2));

    %Yd = [Yd; out(1:2).'];
    Pm = RealMarkerPos(out,Ctr,cam);
    
    Q0 = IK(Pm,shp,Q0);
    %data = [data;Q0.'];
    
    Pix = FK2Pix(Q0,Ctr,cam,shp);
    
    set(h3,'XData',Pix(:,1),'YData',Pix(:,2));
    
    %p = Controller(brd,shp,Q0,G);
    %p = [20, -20]*smoothstep(0.5*brd.t)*cos(0.85*brd.t)^2 * (sign(cos(0.75*brd.t)));
    %p = [20*cos(brd.t),-10*sin(2*brd.t)];
    p = [20,-10];
    
    if brd.t < 5
       %gif('movie_exp.gif','frame',gcf,'nodither');
       %p = [20,-10]*smoothstep(2*brd.t);
       p = lerp([0,0],[20,-10],smoothstep(2*brd.t));
       Pd = [845.6397  325.8011];
    else
       p = lerp([20,-10],[-10,30],smoothstep(0.75*brd.t));
       Pd = [595.3079  319.8983];
       %gif;
    end
    
    if brd.t <= (1/10)
       gif('movie_exp.gif','frame',gcf,'nodither');
    else
       gif;
    end
    
    Pix(end,:)
    
    delete(h4);
    h4 = plot(Pd(1)+x,Pd(2)+y,'m-','LineW',3);
%     
    
    % send data 
    brd.tcpSendData(pdac([p(1), p(2)]));
    
    fprintf('time = %2.3f\n',brd.t);
    drawnow;
end

brd.tcpSendData([0.5, 0.5]);
brd.disconnect();

%% controller
function pd = Controller(brd,shp,Q_,G)
t = brd.t;

% k1 = 1e-4;
% k2 = 1e-2;
% lam1 = 5;
% lam2 = 1e-4;
% Kp = diag([k1,k1,k1,k2,k2,k2]);
% K = 8000*diag([1,1]);
% % 
[g,J] = shp.string(Q_);
% 
gd = SE3(roty(pi/2),[-20,0,-90].');
% 
P   = g(1:3,4,end);
Pd  = gd(1:3,4,end);
% gg = g(:,:,end);
JJ = J(4:6,:,end);
% 
% Xi = smoothstep(t)*logmapSE3(gg\gd);
% Fu = Kp*tmapSE3(Xi)*isomse3(Xi);
% 
% dq = lam1*JJ.'*((JJ*JJ.' + lam2*eye(6))\Fu);
% pd = smoothstep(t)*((G.'*G)/G.')*(K*dq)
%pd = smoothstep(t)*(20*([sin(0.75*pi*t),cos(0.75*pi*t)]) + 10);
pd = G*smoothstep(0.2*t)*1e-3*JJ.'*(Pd-P);
end

function Ctr = CenterMarker(cam)
[x, y] = circlepoints(5);
C = mean(cam.Center,1);
C(1) = C(1) + 12;
C(2) = C(2) + 100;

hold on;
plot(x+C(1), y+C(2),'m-','LineWidth',3);
Ctr = C;
end

function [shp] = SetupModel()
L = 86;  % intrinsic length in mm
M = 1;   % number of modes
N = 11;  % grid on SR
x = linspace(0,L,N).';
Y = zeros(N,M);

for ii = 1:M
    Y(:,ii) = pcc(x/L,ii,M);       % piece-wise constant
end

shp = Shapes(Y,[0,M,0,M,0,0],'L0',L,'g0',SE3(roty(pi/2),zeros(3,1)));
end

function Pix = FK2Pix(q,Ctr,cam,shp)
p = shp.FK(q); P = [p(:,1),-p(:,3)];

Pix = P./cam.P2X(Ctr) + Ctr; 
%hold on; plot(Pix(:,1),Pix(:,2),'o-',...
%    'LineW',2,'Color',col(1));
%P(end,:)
end

function Pix = Point2Pix(P,Ctr,cam)
Pix = P./cam.P2X(Ctr) + Ctr; 
end

function Pm = RealMarkerPos(Out,Ctr,cam)
Pm = (Out(1:2).' - Ctr).*cam.P2X(Ctr);
Pm(:,2) = -Pm(:,2); 
end

function [Q_,p] = IK(Pm,shp,Q_)

Rd = roty(pi/2);
pd = [Pm(1);0;Pm(2)];
gd = SE3(Rd,pd);
e = 10;

k1   = 1e-4;
k2   = 1e-2;
lam1 = 5;
lam2 = 1e-4;
Kp   = diag([k1,k1,k1,k2,k2,k2]);

while norm(e) > 1
    [g,J] = shp.string(Q_);
    
    p  = g(1:3,4,end);
    gg = g(:,:,end);
    JJ = J(:,:,end);
    
    Xi = logmapSE3(gg\gd);
    Fu = Kp*tmapSE3(Xi)*isomse3(Xi);
    
    dq = lam1*JJ.'*((JJ*JJ.' + lam2*eye(6))\Fu);
    e = (p - pd);
   
    Q_  = Q_ + dq;
end

end

