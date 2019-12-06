function x = VersionNumber(Request)

switch(Request)
    case('interface.lib');  x = SorotokiVersion;
    case('preview.lib');    x = SorotokiVersion;
    case('blender.lib');    x = SorotokiVersion;
    case('mesher.lib');     x = SorotokiVersion;
    case('fem.lib');        x = SorotokiVersion;
    case('topology.lib');   x = SorotokiVersion;
    case('former.lib');     x = SorotokiVersion;
    otherwise;              x = SorotokiVersion;
end

end