clr;
%% build Vision class
figure(101);
cam = Vision(1);

%%
for ii = 1:150
    I = snapshot(cam.Cam);
    
    
    disp(ii);
    if ii == 1,
        h = imshow(I);
        gif('testvideo.gif','frame',gcf,'nodither');
        drawnow;
    else
        set(h,'CData',I);
        drawnow limitrate;
       gif; 
    end
end