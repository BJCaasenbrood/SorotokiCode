function color  = col(n)

Colors = pallateV0;

if nargin < 1
    color = vertcat(Colors{:});
else
    color = hex2rgb(Colors{n},1);
end

end

