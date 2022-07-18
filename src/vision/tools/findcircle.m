function [out,im,im0] = findcircle(im0,R,varargin)
%toc

% imshow(im0, 'Border', 'tight'); 
% drawnow;
% f=getframe; 
% im0 = f.cdata;

if size(im0,3) > 1
    im = rgb2gray(im0);
else
    im = im0;
end

%im(200:450,400:900) = 0;
%imshow(im0);
%
e = edge(im,'canny');
h = circle_hough(e, R, 'same', 'normalise');

% get peaks
out = circle_houghpeaks(h, R, ...
                           'nhoodxy', 5,...
                           'nhoodr', 9,...
                           'npeaks', 1);
                       
%[centers, radii, ~] = imfindcircles(im0,R);
% out(1) = centers(1);
% out(2) = centers(2);
% out(3) = radii;
%out = circle_houghpeaks(h,R,'npeaks', Num);
%out = houghcircles(h, R, R);

% figure(101);
% imshow(im0);
if ~isempty(varargin)
Color = barney(4);
for peak = out
     [x, y] = circlepoints(peak(3)); hold on;
     [x2, y2] = circlepoints(peak(3)*0.2); hold on;
     plot(x+peak(1), y+peak(2),'Color',Color(2,:),'LineWidth',6);
     plot(x2+peak(1), y2+peak(2),'Color',Color(2,:),'LineWidth',6);
end
                      
end

end

