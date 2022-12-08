classdef Vision
     
    properties (Access = public)
       Cam;
       Pipe;            % realSense pipeline
       Colorize;        % realSense color
       AR;
       Center;   
       P2X;
       lastFrame;
       Frame;
       output;
       UpdateGain;
       WorldLimits;
    end
    
    properties (Access = private)
       autoConnect; 
       ColorFilter;
       Marker;
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
       %obj.Cam = webcam(1);
       obj.Pipe = realsense.pipeline();
       obj.Colorize = realsense.colorizer();
       cfg = realsense.config();
       %cfg.enable_all_streams();
       cfg.enable_stream(realsense.stream.color);

       profile = obj.Pipe.start(cfg);
       dev = profile.get_device();
       sens = dev.first('depth_sensor');
       sens.set_option(realsense.option.emitter_enabled,0);
    end
   
    obj.NMarker = 1;
    obj.RMarker = 22;
    obj.UpdateGain = 0.3;
    
%     for ii = 2:2:length(varargin)
%         obj.(varargin{ii}) = varargin{ii+1};
%     end

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
%---------------------------------------------------------------------- set
function Vision = setWorld(Vision,xW,yW)

if isempty(Vision.lastFrame)
    Vision = Vision.getframe;
end

Vision.WorldLimits = imref2d(size(Vision.lastFrame),xW,yW);
end
%--------------------------------------------------------- evalution of SDF
function [Vision,img] = getframe(Vision)
    
if isempty(Vision.Cam)
    fs = Vision.Pipe.wait_for_frames();
    %ir = fs.get_infrared_frame();
    %ir = fs.get_fisheye_frame(1);
    ir = fs.get_color_frame();
    %img = reshape(ir.get_data(),848,480).';
    %img = ir.get_data;
    I = ir.get_data();
    img(:,:,1) = reshape(I(1:3:end),640,480).';
    img(:,:,2) = reshape(I(2:3:end),640,480).';
    img(:,:,3) = reshape(I(3:3:end),640,480).';

    %ir_get.data
else
    img = snapshot(Vision.Cam);
end

if ~isempty(Vision.Frame)
    B = Vision.Frame;
    img = img(B(3):B(4),B(1):B(2),:);
end

Vision.lastFrame = img;

end
%--------------------------------------------------------- evalution of SDF
function imgf = getframefiltered(Vision)
    [~,img] = getframe(Vision);
    if ~isempty(Vision.ColorFilter)
        imgf = colorfilter(img,Vision.ColorFilter);
    else
        imgf = img;
    end
end

%--------------------------------------------------------- evalution of SDF
function h = show(Vision)

if isempty(Vision.lastFrame)
    %fs = Vision.Pipe.wait_for_frames();
    %ir = fs.get_infrared_frame(1);
    %img = reshape(ir.get_data(),848,480).';
    Vision = Vision.getframe;
    img = Vision.lastFrame;
else
    img = Vision.lastFrame;
end

figure(101); 
clf;
%if isempty(Vision.WorldLimits)
%    h = imshow(img,'Border','tight'); 
%else
    h = imshow(img,Vision.WorldLimits); 
%end
drawnow;

end
%--------------------------------------------------------- evalution of SDF
function h = showAR(Vision)

    for ii = 1:numel(Vision.AR)
        hold on;
        Y = Vision.AR{ii}; Y = [Y;Y(1,:)];
        C = Vision.Center(ii,:);
       
        plot(Y(:,1),Y(:,2),'LineWidth',3,'Color',col(ii+3));
        plot(C(:,1),C(2),'o','LineWidth',3,'Color',col(ii+3));
    end

end
%--------------------------------------------------------- evalution of SDF
function preview(Vision)
   preview(Vision.Cam); 
end
%--------------------------------------------------------- evalution of SDF
function [Vision] = detect(Vision,varargin)
CUT = [];    
if isempty(varargin)
    R = Vision.RMarker;
elseif numel(varargin) == 1
    R = varargin{1};
else
    R = varargin{1};
    CUT = varargin{2};
end

if isempty(Vision.lastFrame)
Vision = getframe(Vision);
end
    
im0 = Vision.lastFrame;
%if ~isempty(Vision.Frame)
%   im0 = im0(CUT(1):CUT(2),CUT(3):CUT(4),:); 
%end
%toc

% imshow(im0, 'Border', 'tight'); 
% drawnow;
% f=getframe; 
% im0 = f.cdata;

if ~isempty(Vision.ColorFilter)
    im0 = colorfilter(im0,Vision.ColorFilter);
end

if numel(size(im0)) > 2
    im = rgb2gray(im0);
else
    im = im0;
end

e = edge(im,'canny');
h = circle_hough(e, R,'same','normalise');

% get peaks
out = circle_houghpeaks(h, R, ...
                           'nhoodxy', Vision.HoodXY,...
                           'nhoodr', Vision.HoodR,...
                           'npeaks', Vision.NMarker);

R   = mean(out(3,:));
out = out(1:2,:);

if ~isempty(Vision.WorldLimits)
    [out(1),out(2)] = Vision.WorldLimits.intrinsicToWorld(out(1),out(2));
end

% imshow(im0);
% Color = barney(4);
% 
% for peak = out
%      [x, y] = circlepoints(R); hold on;
%      
%      if ~isempty(Vision.Marker)
%          delete(Vision.Marker);
%      end
% 
%      X = x+peak(1);
%      Y = y+peak(2);
% 
%      if ~isempty(Vision.WorldLimits)
%         [X,Y] = Vision.WorldLimits.intrinsicToWorld(X,Y);
%      end
% 
%      Vision.Marker = plot(X,Y,'Color',Color(2,:),...
%          'LineWidth',3);
%      drawnow();
% end                

if ~isempty(Vision.output)
    Vision.output = lerp(out,Vision.output,Vision.UpdateGain);
else
    Vision.output = out;
end
    
end
%--------------------------------------------------------- evalution of SDF
function [X,Y] = getMarker(Vision,R)
 [x, y] = circlepoints(R); 

 if ~isempty(Vision.Marker)
     delete(Vision.Marker);
 end

 peak = Vision.output;
 X = x+peak(1);
 Y = y+peak(2);

 if ~isempty(Vision.WorldLimits)
     [X,Y] = Vision.WorldLimits.intrinsicToWorld(X,Y);
 end

end

%--------------------------------------------------------- evalution of SDF
function Vision = detectAR(Vision,varargin)
% % get tight image
pathCD = cdsoro;
pathFile = '/src/vision/tools/python/';

if isempty(varargin) 
    %img = snapshot(Vision.Cam);
    if isempty(Vision.lastFrame)
        Vision = Vision.getframe();
    end
    img = Vision.lastFrame;
    imwrite(img,[pathCD,pathFile,'tmp.png']);
end

%system(['cd ',pathCD,pathFile,' && python3 findMarkers.py']);
pyrunfile([pathCD,pathFile,'findMarkers.py']);

in  = importdata(['~',pathCD,pathFile,'readme.txt']);
nAr = numel(in)/8;

cout(['* Detected #AR markers = ', num2str(nAr), '\n']);

Vision.AR = {};
Vision.P2X = [];    
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

