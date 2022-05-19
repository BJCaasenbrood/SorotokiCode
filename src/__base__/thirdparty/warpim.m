function handle = warpim(varargin)
%WARP Display image as texture-mapped surface.
%   WARP(X,MAP) displays the indexed image X with colormap MAP as
%   a texture map on a simple rectangular surface.
%
%   WARP(I,N) displays the intensity image I with gray scale
%   colormap of length N as a texture map on a simple rectangular
%   surface.
%
%   WARP(BW) displays the binary image BW as a texture map on a
%   simple rectangular surface.
%
%   WARP(RGB) displays the RGB image in the array RGB as a
%   texture map on a simple rectangular surface.
%
%   WARP(z,...) displays the image on the surface z.
%
%   WARP(x,y,z,...) displays the image on the surface (x,y,z).
%
%   H = WARP(...) returns a handle to the texture mapped
%   surface.
%
%   Class Support
%   -------------
%   The input image can be of class logical, uint8, uint16, or double.
%
%   Remarks
%   -------
%   Texture-mapped surfaces generally render more slowly than
%   images.
%
%
%   See also IMSHOW, IMAGE, IMAGESC, SURF.

%   Copyright 1993-2016 The MathWorks, Inc.

[x,y,z,cdata,cdatamapping,clim,map,likeimage] = parse_inputs(varargin{:});

axHandle = newplot;
h = surface(x,y,z,cdata,'EdgeColor','none','FaceColor','texturemap', ...
        'CDataMapping',cdatamapping);
if (~isempty(clim))
    set(axHandle, 'CLim', clim);
end
if (~isempty(map))
    axHandle.ColorSpace.Colormap = map;
end

if likeimage && ~ishold
    %view(2)
    %axis([min(x(:)) max(x(:)) min(y(:)) max(y(:))])
else
    %view(3)
end

if nargout, handle = h; end


%-----------------------------------------------------------
% Subfunction PARSE_INPUTS
%-----------------------------------------------------------
function [x,y,z,cdata,cdatamapping,clim,map,likeimage] = ...
        parse_inputs(varargin)

x = [];
y = [];
z = [];
map = [];
cdatamapping = 'direct';
clim = [];
likeimage = 0;
if (get(0,'ScreenDepth') > 16)
    defGrayMapLength = 256;
else
    defGrayMapLength = 64;
end

switch nargin
case 0
    error(message('images:warp:notEnoughInputs'))
    
case 1
    % warp(I)
    % warp(RGB)
    likeimage = 1;

    if ((ndims(varargin{1}) == 3) && (size(varargin{1},3) == 3))
        % warp(RGB)
        cdata = varargin{1};
        if (~isa(cdata,'double'))
            cdata = im2double(cdata);
        end
        
    else
        % warp(I)
        cdata = varargin{1};
        cdatamapping = 'scaled';
        clim = [0 1];
        if (~isa(cdata,'double'))
            cdata = im2double(cdata);
        end
        map = gray(defGrayMapLength);
    end
    
case 2
    % warp(X,map)
    % warp(I,N)
    % warp(z,I)
    % warp(z,RGB)
    % warp(I,[a b])
    
    if ((ndims(varargin{2}) == 3) && (size(varargin{2},3) == 3))
        % warp(z,RGB)
        z = varargin{1};
        cdata = varargin{2};
        if (~isa(cdata,'double'))
            cdata = im2double(cdata);
        end
        
    elseif (numel(varargin{2}) == 1)
        % warp(I,N)
        cdata = varargin{1};
        map = gray(varargin{2});
        cdatamapping = 'scaled';
        clim = [0 1];
        if (~isa(cdata,'double'))
            cdata = im2double(cdata);
        end
        likeimage = 1;
        
    elseif (isequal(size(varargin{2}), [1 2]))
        % warp(I,[a b]) 
        cdata = varargin{1};
        cdatamapping = 'scaled';
        clim = varargin{2};
        map = gray(defGrayMapLength);
        
        if isa(cdata,'uint8')
            cdata = im2double(cdata);
            clim = clim/255.0;
        elseif isa(cdata,'uint16')
            cdata = im2double(cdata);
            clim = clim/65535.0;
        elseif islogical(cdata)
            cdata = im2double(cdata);
        end
            
        likeimage = 1;
        
    elseif (size(varargin{2},2) == 3)
        % warp(X,map)
        cdata = varargin{1};
        map = varargin{2};
        cdatamapping = 'direct';
        if ~isa(cdata,'double')
            cdata = im2double(cdata, 'indexed');
        end
        likeimage = 1;
        
    else
        % warp(z,I)
        z = varargin{1};
        cdata = varargin{2};
        cdatamapping = 'scaled';
        clim = [0 1];
        if ~isa(cdata,'double')
            cdata = im2double(cdata, 'indexed');
        end
        map = gray(defGrayMapLength);
        
    end
    
case 3
    % warp(z,X,map)
    % warp(z,I,N)
    % warp(z,I,[a b])
    
    if (numel(varargin{3}) == 1)
        % warp(z,I,N)
        z = varargin{1};
        cdata = varargin{2};
        map = gray(varargin{3});
        cdatamapping = 'scaled';
        clim = [0 1];
        if ~isa(cdata,'double')
            cdata = im2double(cdata);
        end
        
    elseif (isequal(size(varargin{3}), [1 2]))
        % warp(z,I,[a b])
        z = varargin{1};
        cdata = varargin{2};
        cdatamapping = 'scaled';
        clim = varargin{3};
        map = gray(defGrayMapLength);
        if isa(cdata,'uint8')
            cdata = im2double(cdata);
            clim = clim/255.0;
        elseif isa(cdata,'uint16')
            cdata = im2double(cdata);
            clim = clim/65535.0;
        elseif islogical(cdata)
            cdata = im2double(cdata);
        end
        
    elseif (size(varargin{3},2) == 3)
        % warp(z,X,map)
        z = varargin{1};
        cdata = varargin{2};
        map = varargin{3};
        cdatamapping = 'direct';
        if ~isa(cdata,'double')
            cdata = im2double(cdata, 'indexed');
        end
        
    else
        error(message('images:warp:invalidInputs'))
        
    end
    
case 4
    % warp(x,y,z,I)
    % warp(x,y,z,RGB)
    
    if ((ndims(varargin{4}) == 3) && (size(varargin{4},3) == 3))
        % warp(x,y,z,RGB)
        x = varargin{1};
        y = varargin{2};
        z = varargin{3};
        cdata = varargin{4};
        if ~isa(cdata,'double')
            cdata = im2double(cdata);
        end
        
    else
        % warp(x,y,z,I)
        x = varargin{1};
        y = varargin{2};
        z = varargin{3};
        cdata = varargin{4};
        cdatamapping = 'scaled';
        clim = [0 1];
        map = gray(defGrayMapLength);
        if ~isa(cdata,'double')
            cdata = im2double(cdata);
        end
        
    end
    
case 5
    % warp(x,y,z,X,map)
    % warp(x,y,z,I,N)
    % warp(x,y,z,I,[a b])
    
    if (numel(varargin{5}) == 1)
        % warp(x,y,z,I,N)
        x = varargin{1};
        y = varargin{2};
        z = varargin{3};
        cdata = varargin{4};
        map = gray(varargin{5});
        cdatamapping = 'scaled';
        clim = [0 1];
        if ~isa(cdata,'double')
            cdata = im2double(cdata);
        end
        
    elseif (isequal(size(varargin{5}), [1 2]))
        % warp(x,y,z,I,[a b])
        x = varargin{1};
        y = varargin{2};
        z = varargin{3};
        cdata = varargin{4};
        cdatamapping = 'scaled';
        clim = varargin{5};
        map = gray(defGrayMapLength);
        if isa(cdata,'uint8')
            cdata = im2double(cdata);
            clim = clim/255.0;
        elseif isa(cdata,'uint16')
            cdata = im2double(cdata);
            clim = clim/65535.0;
        elseif islogical(cdata)
            cdata = im2double(cdata);
        end
        
    elseif (size(varargin{5},2) == 3)
        % warp(x,y,z,X,map)
        x = varargin{1};
        y = varargin{2};
        z = varargin{3};
        cdata = varargin{4};
        map = varargin{5};
        cdatamapping = 'direct';
        if ~isa(cdata,'double')
            cdata = im2double(cdata, 'indexed');
        end
        
    else
        error(message('images:warp:invalidInputs'))
        
    end
    
otherwise
    error(message('images:warp:tooManyInputs'))
    
end

siz = size(cdata);
M = siz(1);
N = siz(2);

if (isempty(z))
    % The surface displays most quickly when we use
    % a simple 2-by-2 z matrix, but that uses up a
    % large amount of printer memory when printed.
    % In IPT v1, the z matrix was the same size as
    % the image; this solution took a long time to
    % display.  The factor of 4 below is a compromise.
    % -sle, September 1996
    p = max(floor(min(size(cdata))/4),2);
    z = zeros(p);
    x = linspace(1,N,p);
    y = linspace(1,M,p);
end

if (isempty(x))
    [x,y] = meshgrid(1:size(z,2), 1:size(z,1));
end

if ((length(x) == 2) && (length(y) == 2))
    [x,y] = meshgrid(linspace(x(1),x(2),size(z,2)), ...
            linspace(y(1),y(2),size(z,1)));
end
