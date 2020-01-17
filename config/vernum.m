function x = vernum(Request)

switch(Request)
    case('base.lib');       x = '3.0';
    case('gmdl.lib');       x = '3.0';
    case('mesh.lib');       x = '3.0';
    case('fem.lib');        x = '3.0';
    case('mdl.lib');        x = '3.0';
    otherwise;              x = 0;
end

end