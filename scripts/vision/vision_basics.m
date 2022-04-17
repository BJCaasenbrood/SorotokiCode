clr;
%% build Vision class
cam = Vision();

%% detect cirlces
out = cam.detect();

%% detect AR markers
h = cam.show();
tic;
cam = cam.detectAR();
toc;

for ii = 1:length(cam.AR) 
   Nds = cam.AR{ii};
   hold on; plot(Nds(:,1),Nds(:,2),...
       '-o','LineW',5,'Color',col(4));
   plot(cam.Center(ii,1),cam.Center(ii,2),'Color',col(4),...
       'LineW',3);
   
   [x, y] = circlepoints(out(3));
   plot(x+out(1), y+out(2),...
        'g-','LineWidth',3);
   plot(out(1),out(2),'go','LineW',3); 
   plot(cam.Center(ii,1),cam.Center(ii,2),'bo','LineW',3); 
end

[x, y] = circlepoints(5);
C = mean(cam.Center,1);
C(1) = C(1) + 7;
C(2) = C(2) + 100;
plot(x+C(1), y+C(2),...
        'm-','LineWidth',3);