classdef Vision
     
    properties (Access = public)
       Cam;
       Pipe;    %real sense pipeline
       Colorize;
       AR;
       Center;   
       P2X;
    end
    
    properties (Access = private)
       autoConnect; 
       NMarker;
       RMarker;
       HoodXY = 5;
       HoodR  = 9;
    end
%--------------------------------------------------------------------------    
%--------------------------------------------------------------------------       
methods (Access = public)
%--------------------------------------------------------------------------       
function obj = Vision(varargin)
    
    list = webcamlist;
    
    if isempty(varargin)
       Id = action(list);
       obj.Cam = webcam(list{Id});
    elseif isfloat(varargin{1})
       Id = varargin{1};
       varargin{1} = [];
       obj.Cam = webcam(list{Id});
    elseif strcmp('realsense',varargin{1})
       obj.Cam = webcam(1);
       obj.Pipe = realsense.pipeline();
       obj.Colorize = realsense.colorizer();
    end
    
%         for ii = 1:2:length(varargin)
%             obj.(varargin{ii}) = varargin{ii+1};
%         end
%     
    obj.NMarker = 1;
    obj.RMarker = 22;
    
end
%---------------------------------------------------------------------- get     
function varargout = get(Vision,varargin)
    if nargin > 1
        varargout{nargin-1,1} = [];
        for ii = 1:length(varargin)
            varargout{ii,1} = Vision.(varargin{ii});
        end
    else
        varargout = Vision.(varargin);
    end
end
%---------------------------------------------------------------------- set
function Vision = set(Vision,varargin)
    for ii = 1:2:length(varargin)
        Vision.(varargin{ii}) = varargin{ii+1};
    end
end
%--------------------------------------------------------- evalution of SDF
function h = show(Vision)
    
img = snapshot(Vision.Cam);
figure(101); 
clf;
h = imshow(img,'Border','tight'); 
drawnow;
% f    = getframe; 
% imss = f.cdata;

end
%--------------------------------------------------------- evalution of SDF
function preview(Vision)
   preview(Vision.Cam); 
end
%--------------------------------------------------------- evalution of SDF
function [out,im,im0] = detect(Vision,varargin)
CUT = [];    
if isempty(varargin)
    R = Vision.RMarker;
else
    R = varargin{1};
    CUT = varargin{2};
end



Num = Vision.NMarker;
%tic
im0 = snapshot(Vision.Cam);
if ~isempty(CUT)
   im0 = im0(CUT(1):CUT(2),CUT(3):CUT(4),:); 
end
%toc

% imshow(im0, 'Border', 'tight'); 
% drawnow;
% f=getframe; 
% im0 = f.cdata;

im = rgb2gray(im0);

%im(200:450,400:900) = 0;
%imshow(im0);
%
e = edge(im,'canny');
h = circle_hough(e, R, 'same', 'normalise');

% get peaks
out = circle_houghpeaks(h, R, ...
                           'nhoodxy', Vision.HoodXY,...
                           'nhoodr', Vision.HoodR,...
                           'npeaks', Num);
                       
%[centers, radii, ~] = imfindcircles(im0,R);
% out(1) = centers(1);
% out(2) = centers(2);
% out(3) = radii;
%out = circle_houghpeaks(h,R,'npeaks', Num);
%out = houghcircles(h, R, R);

figure(101);
imshow(im0);
Color = barney(4);
for peak = out
     [x, y] = circlepoints(peak(3)); hold on;
     [x2, y2] = circlepoints(peak(3)*0.2); hold on;
     plot(x+peak(1), y+peak(2),'Color',Color(2,:),'LineWidth',3);
     plot(x2+peak(1), y2+peak(2),'Color',Color(2,:),'LineWidth',3);
end
                      
    
end
%--------------------------------------------------------- evalution of SDF
function Vision = detectAR(Vision,varargin)
% % get tight image
pathCD = cdsoro;
pathFile = '/src/vision/tools/python/';

if isempty(varargin) 
    img = snapshot(Vision.Cam);
    imwrite(img,[pathCD,pathFile,'tmp.png']);
end

system(['cd ',pathCD,pathFile,' && python3 findMarkers.py']);

in  = importdata([pathCD,pathFile,'readme.txt']);
nAr = numel(in)/8;

cout(['* Detected #AR markers = ',num2str(nAr), '\n']);

Vision.AR = {};
for ii = 1:nAr
    Vision.AR{ii,1}     = in(4*ii-3:4*ii,:);
    Vision.Center(ii,:) = mean(in(4*ii-3:4*ii,:),1);
    
    dX = sqrt(sum(diff(Vision.AR{ii,1}(1:3,:)).^2,2));
    Vision.P2X(ii,:)    = 30./dX; 
end
    
C = Vision.Center;
G = Vision.P2X;

Vision.P2X = @(x) [lerp(G(1,1),G(1,2),invlerp(C(1,1),C(1,2),x(1))),...
                     lerp(G(2,1),G(2,2),invlerp(C(2,1),C(2,2),x(2)))];
end
end

%--------------------------------------------------------- evalution of SDF
methods (Access = private)

end
end

