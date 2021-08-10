function color  = col(n,per)

Colors = pallateV0;

if nargin < 1
    C = vertcat(Colors{:});
    color = [];
    for ii = 1:length(C)
        color = [color;hex2rgb(Colors{ii},1)];
    end
else
    color = hex2rgb(Colors{n},1);
end

if nargin > 1
    ColorNew = SrgbCompanding(color,0);
    ColorNull = SrgbCompanding([1,1,1],0);
    color = AffineColorMix(ColorNew,ColorNull,[0,0,0],[(per),...
        (1-(per)),0]); 
    color = SrgbCompanding(color,1);
end

end

function Colors = pallateV0
Colors{1}     = '#0058A8'; % blue-soft 0053D0 0161E8 4169E1
Colors{end+1} = '#CA1E2C'; % red-dark :was E82C0C
Colors{end+1} = '#34A73F'; % green-dark
Colors{end+1} = '#EE7023'; % yellow FF9632 FF901F
Colors{end+1} = '#7E2583'; % purple-dark
Colors{end+1} = '#F9AB15';
Colors{end+1} = '#3E2E86'; % pink-dark
Colors{end+1} = '#B6D124'; % pink-dark F4A460 F222FF
Colors{end+1} = '#5F2A84'; % blue
Colors{end+1} = '#E64928';
Colors{end+1} = '#FFF100'; % pink-dark
Colors{end+1} = '#009D7E'; % pink-dark F4A460 F222FF
end

function Colors = pallateV1
Colors{1}     = '#0053D0'; %  blue-soft 0053D0 0161E8 4169E1
Colors{end+1} = '#E62014'; % red-dark :was E82C0C
Colors{end+1} = '#5CB52F'; % green-dark
Colors{end+1} = '#FF8B17'; % yellow FF9632 FF901F
Colors{end+1} = '#9B00E8'; % purple-dark
Colors{end+1} = '#0265FD';
Colors{end+1} = '#FF2975'; % pink-dark
Colors{end+1} = '#F4A460'; % pink-dark F4A460 F222FF
Colors{end+1} = '#290CFF'; % blue
end

function  Colors = pallateV2
Colors{1}     = '#1b61a4'; % blue-soft
Colors{end+1} = '#f51e3d'; % red-soft
Colors{end+1} = '#27ac4d'; % green-soft
Colors{end+1} = '#f1621e'; % orange-soft
Colors{end+1} = '#ec237d'; % purple-dark
Colors{end+1} = '#eee70b'; % bright-yellow
Colors{end+1} = '#df0a36'; % dark-rouge-dark
Colors{end+1} = '#71ba8d'; % ligth-green
Colors{end+1} = '#290CFF'; % blue
end

function rgb = hex2rgb(hex,range)
% * * * * * * * * * * * * * * * * * * * * 
% Chad A. Greene, April 2014
% * * * * * * * * * * * * * * * * * * * * 
% Input checks:
assert(nargin>0&nargin<3,'hex2rgb function must have one or two inputs.') 
if nargin==2
    assert(isscalar(range)==1,'Range must be a scalar, either "1" to scale from 0 to 1 or "256" to scale from 0 to 255.')
end
% Tweak inputs if necessary: 
if iscell(hex)
    assert(isvector(hex)==1,'Unexpected dimensions of input hex values.')

    if isrow(hex)
        hex = hex'; 
    end
    
    % If input is cell, convert to matrix: 
    hex = cell2mat(hex);
end
if strcmpi(hex(1,1),'#')
    hex(:,1) = [];
end
if nargin == 1
    range = 1; 
end
% Convert from hex to rgb: 
switch range
    case 1
        rgb = reshape(sscanf(hex.','%2x'),3,[]).'/255;
    case {255,256}
        rgb = reshape(sscanf(hex.','%2x'),3,[]).';
    
    otherwise
        error('Range must be either "1" to scale from 0 to 1 or "256" to scale from 0 to 255.')
end

end

%------------------------------------------- CONVERT/INVERT sRGB COMPANDING
function ColorNew = SrgbCompanding(Color,Inverse)
r = Color(1);
g = Color(2);
b = Color(3);

if Inverse == 1
    if (r > 0.04045), r = ((r+0.055)/1.055)^2.4; else, r = r / 12.92; end
    if (g > 0.04045), g = ((g+0.055)/1.055)^2.4; else, g = g / 12.92; end
    if (b > 0.04045), b = ((b+0.055)/1.055)^2.4; else, b = b / 12.92; end
else
    if (r > 0.0031308),r = 1.055*(r)^(1/2.4)-0.055; else, r = r * 12.92;end
    if (g > 0.0031308),g = 1.055*(g)^(1/2.4)-0.055; else, g = g * 12.92;end
    if (b > 0.0031308),b = 1.055*(b)^(1/2.4)-0.055; else, b = b * 12.92;end
end

ColorNew = [r, g, b];
end

%------------------------------------------------ AFFINE COLOR COMBINATIONS
function Color = AffineColorMix(C1,C2,C3,Mix)
a1 = Mix(1);
a2 = Mix(2);
a3 = Mix(3);

Color1 = SrgbCompanding(C1,1);
Color2 = SrgbCompanding(C2,1);
Color3 = SrgbCompanding(C3,1);

r = Color1(1)*a1 + Color2(1)*a2 + Color3(1)*a3;
g = Color1(2)*a1 + Color2(2)*a2 + Color3(2)*a3;
b = Color1(3)*a1 + Color2(3)*a2 + Color3(3)*a3;

Color = [r,g,b];

Color = SrgbCompanding(Color,0);
end

