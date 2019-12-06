function x = vernum(Request)

switch(Request)
    case('base.lib');       x = '2.1';
    case('output.lib');     x = '2.1';
    case('object.lib');     x = '2.1';
    otherwise;              x = 0;
end

end