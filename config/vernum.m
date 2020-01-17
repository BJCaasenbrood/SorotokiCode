function x = vernum(Request)

switch(Request)
    case('base.lib');       x = '2.1';
    case('gmdl.lib');       x = '2.1';
    case('mesh.lib');       x = '2.1';
    case('fem.lib');        x = '2.1';
    case('mdl.lib');        x = '2.1';
    otherwise;              x = 0;
end

end