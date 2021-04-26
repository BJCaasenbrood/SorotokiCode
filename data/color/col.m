function color  = col(n,per)
Colors = pallateV1;
if nargin < 1
    color = vertcat(Colors{:});
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





