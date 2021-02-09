function color  = col(n,per)
Colors = pallateV0;
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

