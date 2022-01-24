function h = shade(varargin)
%SHADE  Filled area linear plot.
%   
%   SHADE should be called using the same syntax as the built-in PLOT.
%
%   SHADE(X,Y) plots vector Y versus vector X, filling the area under the
%   curve. If either is a matrix, this function behaves like PLOT.
%
%   SHADE(Y) plots the columns of Y versus their index.
%   
%   SHADE(X,Y,S) plots Y versus X using the line type, marker symbols and
%   colors as specified by S. For more information on the line specifier S,
%   see PLOT.
%
%   SHADE(X1,Y1,X2,Y2,X3,Y3,...) combines the plots defined by
%   the (X,Y) pairs.
%
%   SHADE(X1,Y1,S1,X2,Y2,S2,X3,Y3,S3,...) combines the plots defined by
%   the (X,Y,S) triples.
%
%   SHADE(AX,...) plots into the axes with handle AX.
%
%   H = SHADE(...) returns a column vector of handles to graphics objects.
%
%   SHADE(...,Name,Value) specifies additional properties of the lines (see
%   help for PLOT) or the filled areas (see below).
%
%   Three additional properties are provided to control the filling:
%
%     - 'FillType' specifies the filling behaviour. The input should be a
%       matrix of size [N,2], such as [A1,B1;A2,B2;...;An,Bn], where Ai and
%       Bi indicate the upper and lower limits, respectively, of each of
%       the N areas to be filled. Each Ai and Bi is an index pointing to
%       one of the lines drawn by PLOT. If, for a particular combination of
%       inputs, PLOT draws M lines, then each Ai and Bi should be a number
%       between 1 and M. In addition, the special cases 0, -1 and -2 are
%       allowed, each representing the x-axis, the bottom of the active
%       axes and the top of the active axes, respectively. 'FillType' may
%       also be specified using a cell array of size [N,2], in which case
%       one may also use the keywords 'axis', 'bottom' and 'top' instead
%       of 0, -1 and -2. By default, the areas between each of the curves
%       and the x-axis are filled, which corresponds to the input matrix
%       [1,0;0,1;2,0;0,2;...;M,0;0,M].
%
%     - 'FillColor' specifies the color of the fillings. This should be a
%       matrix of size [N,3], where each row is the RGB triplet for the
%       corresponding area as specified above. 'FillColor' may also be
%       specified as a cell array of length N, in which case each entry may
%       be either an RGB triplet or any of the color names commonly used in
%       MATLAB. If only one RGB value or color name is provided, all areas
%       are treated equally. If this parameter is not specified, colors are
%       determined by the corresponding lines.
%
%     - 'FillAlpha' specifies the transparency of the fillings. This should
%       be a vector of length N, where each entry specifies the alpha value
%       for each area. If only one alpha value is provided, all areas are
%       treated equally. If this parameter is not specified, an alpha value
%       of 0.3 is used for all areas.
%
%   See also PLOT.

% Copyright (c) 2018 Javier Montalt Tordera.

% accepted params
names = {'FillType','FillColor','FillAlpha'};

% init fill params
fp = cell(1,3);

% extract filling parameters, if present
for n = 1:length(names)
    for i = 1:length(varargin)
        % if found
        if strcmpi(names{n},varargin{i})
            if i+1 > nargin
                error(["Expected an input value after the name '" names{i} "'."]);
            end
            % save filling info
            fp{n} = varargin{i+1};
            % delete from varargin array - otherwise plot will fail as it won't
            % understand the input
            varargin(i:i+1) = [];
            break;
        end
    end
end

% check if an axes object was specified
if isscalar(varargin{1}) && ishandle(varargin{1}(1))
    ax = varargin{1};
else
    ax = gca;
end

% initial hold state
tf = ishold(ax);

% plot lines
ls = plot(varargin{:});
hold(ax,'on');

% provide default filling params
fd = fp{1};
fc = fp{2};
fa = fp{3};

% validate fill type
fd = validatetype(fd,ls);

% number of fillings
nf = size(fd,1);

% validate fill color
fc = validatecolor(fc,fd,ls,nf);

% validate fill alpha
fa = validatealpha(fa,nf);

% array to hold patch objects
ps = gobjects(nf,1);

% for each filling
for i = 1:nf
    
    x = cell(1,2);
    y = cell(1,2);
    
    % get data
    for j = 1:2
        switch fd(i,j)
            case {-2,-1}
                y{j} = ylim;
                y{j} = y{j}(abs(fd(i,j)));
            case 0
                y{j} = 0;
            otherwise
                x{j} = ls(fd(i,j)).XData;
                y{j} = ls(fd(i,j)).YData;
        end
    end
    
    if isequal(x{1},x{2})
        x = x{1};
    elseif isempty(x{1})
        x = x{2};
        y{1} = y{1} * ones(size(x));
    elseif isempty(x{2})
        x = x{1};
        y{2} = y{2} * ones(size(x));
    else
        x = sort([x{1} x{2}]);
        y{1} = interp1(x{1},y{1},x,'linear','extrap');
        y{2} = interp1(x{2},y{2},x,'linear','extrap');
    end
    
    % crossings
    x0 = [x(1) zcross(x,y{1} - y{2}) x(end)];
    y0 = interp1(x,y{1},x0);
    
    % for each zero-crossing
    for j = 0:length(x0)-2

        % index
        idx = x >= x0(j+1) & x <= x0(j+2) & y{1} >= y{2};
        if all(~idx), continue; end

        % polygon corners
        xv = [x0(j+1) x(idx) x0(j+2) fliplr(x(idx))];
        yv = [y0(j+1) y{1}(idx) y0(j+2) fliplr(y{2}(idx))];

        % fill polygon
        ps(i) = fill(ax,xv,yv,fc(i,:),'LineStyle','none','FaceAlpha',fa(i));

    end
end

% release if the hold state was off when called
if tf == 0
    hold(ax,'off');
end

% set output argument if requested
if nargout == 1
    h = [ls;ps];
end

end


function z = zcross(x,y)
% find zero crossings of line Y versus X

% logical index
c = y > 0;

% find point pairs where there is a sign change
d = abs(diff(c));
p1 = find(d == 1);  % before change
p2 = p1 + 1;        % after change

% zero-crossing positions
z = x(p1) + abs(y(p1)) ./ (abs(y(p1)) + abs(y(p2))) .* (x(p2) - x(p1));

end


function fd = validatetype(fd,ls)

% if no filling specified, fill to x-axis
if isempty(fd)
    fd = [(1:length(ls))' zeros(size(ls))];
    fd = [fd;fliplr(fd)];
    fd = fd([1:2:length(fd) 2:2:length(fd)],:);
end

% if the filling was specified in cell form, convert to matrix
if iscell(fd)
    tmp = zeros(size(fd));
    for i = 1:numel(fd)
        if ischar(fd{i})
            tmp(i) = str2ind(validatestring(fd{i},{'axis','bottom','top'},'shade','FillType'));
        else
            validateattributes(fd{i},{'numeric'},{'scalar'},'shade','FillType');
            tmp(i) = fd{i};
        end
    end
    fd = tmp;
end

validateattributes(fd,{'numeric'},{'integer','size',[nan 2],'>=',-2,'<=',length(ls)},'shade','FillType');

end


function fc = validatecolor(fc,fd,ls,nf)

% if no color specified, get it from plot lines
if isempty(fc)
    fc = zeros(nf,3);
    for i = 1:nf
        if fd(i,1) <= 0
            fc(i,:) = ls(fd(i,2)).Color;
        elseif fd(i,2) <= 0
            fc(i,:) = ls(fd(i,1)).Color;
        else
            fc(i,:) = mean([ls(fd(i,1)).Color;ls(fd(i,2)).Color],1);
        end
    end
end

% if length 1, repeat
if ischar(fc)
    fc = {fc};
    fc = repmat(fc,nf,1);
elseif (iscell(fc) && numel(fc) == 1) || (~iscell(fc) && size(fc,1) == 1)
    fc = repmat(fc,nf,1);
end

% if color is specified in cell form, convert to matrix
if iscell(fc)
    validateattributes(fc,{'cell'},{'vector','numel',nf},'shade','FillColor');
    tmp = zeros(numel(fc),3);
    for i = 1:numel(fc)
        if ischar(fc{i})
            tmp(i,:) = str2rgb(validatestring(fc{i},{'y','m','c','r','g','b','w','k','yellow','magenta','cyan','red','green','blue','white','black'},'shade','FillColor'));
        else
            validateattributes(fc{i},{'numeric'},{'vector','numel',3},'shade','FillColor');
            if iscolumn(fc{i})
                fc{i} = fc{i}';
            end
            tmp(i,:) = fc{i};
        end
    end
    fc = tmp;
end     

validateattributes(fc,{'numeric'},{'real','size',[nf 3],'>=',0,'<=',1},'shade','FillColor');

end


function fa = validatealpha(fa,nf)

% if no alpha specified, choose a value of 0.2
if isempty(fa)
    fa = 0.3 * ones(nf,1);
end

% if length 1, repeat
if length(fa) == 1
    fa = repmat(fa,nf,1);
end

if isrow(fa)
    fa = fa';
end

validateattributes(fa,{'numeric'},{'real','vector','numel',nf,'>=',0,'<=',1},'shade','FillAlpha');

end


function n = str2ind(s)
% convert string to index

switch s
    case 'axis'
        n = 0;
    case 'bottom'
        n = -1;
    case 'top'
        n = -2;
end

end


function n = str2rgb(s)
% convert string to RBG triplet

switch s
    case {'y','yellow'}
        n = [1 1 0];
    case {'m','magenta'}
        n = [1 0 1];
    case {'c','cyan'}
        n = [0 1 1];
    case {'r','red'}
        n = [1 0 0];
    case {'g','green'}
        n = [0 1 0];
    case {'b','blue'}
        n = [0 0 1];
    case {'w','white'}
        n = [1 1 1];
    case {'k','black'}
        n = [0 0 0];
end


end
