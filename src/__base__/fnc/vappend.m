function y = vappend(vec,x)
[s1,s2] = size(vec);

if (s2 == 1) && (length(x) == 1), y = zeros(s1+1,1); flag = 1;
elseif (s1 == 1) && (length(x) == 1), y = zeros(1,s2+1); flag = 2;
elseif (s1 == length(x)) && (s2 == length(x)), y = zeros(s1+1,s2); flag = 3;
elseif s1 == length(x), y = zeros(s1,s2+1); flag = 4;
elseif s2 == length(x), y = zeros(s1+1,s2); flag = 5;
else, error('Dimensions are incorrect');
end

switch(flag)
    case(1), y = [vec;x];
    case(2), y = [vec,x];
    case(3), y = [vec; x(:)'];
    case(4), y = [vec,x(:)];
    case(5), y = [vec; x(:)'];
end

end

