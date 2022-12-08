clr;
%h = animatedline('MaximumNumPoints',100);
cam = Vision('realsense');
cam = cam.set('HoodXY',21,'HoodR',51,...
    'Frame',[190,553,137,409],...
    'ColorFilter',[360,50],...
    'UpdateGain',0.4);

h = []; h2 = [];

Y = [];
while true
    tic;
    
    delete(h2);
    %addpoints(h,out(2),out(1));
    cam = cam.getframe();
    if isempty(h)
        h = imshow(cam.lastFrame);
    else
        h.CData = cam.lastFrame;
    end

    cam = cam.detect(9);
   
    Y = [Y;cam.output(1:2).'];

    drawnow;
    t = toc;
end

                