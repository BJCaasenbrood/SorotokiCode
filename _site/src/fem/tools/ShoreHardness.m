function E = ShoreHardness(value,SCALE)
%
switch(SCALE)
    case('A'); E = ShoreA(value);             
    case('D'); E = ShoreD(value);             
    otherwise; E = ShoreA(value);             
end 

end

%https://en.wikipedia.org/wiki/Shore_durometer
% x = ASTM D2240 type A hardness.
function E = ShoreA(S)
if (S<20) || (S>80)
    warning('Shore hardness must be 20 < Sd < 80') ;
end
E = fzero(@(x) log(x)/log(10) - 0.0235*S - 0.6403,150);
end

function E = ShoreD(S)
if (S<30) || (S>85)
warning('Shore hardness must be 30 < Sd < 85') ;
end

E = fzero(@(x) log(x)/log(10) - 0.0235*S - 0.6403,150);
end
