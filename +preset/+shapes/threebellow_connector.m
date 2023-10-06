function [obj1, obj2, varargout] = threebellow_connector(varargin)

    p = inputParser;
    addOptional(p,'n',10);
    addOptional(p,'m',1);
    addOptional(p,'phi',0);
    addOptional(p,'q0',1e-2);
    parse(p,varargin{:});

    currentDir = fileparts(mfilename('fullpath'));
    stlPath1 =fullfile(currentDir, '..', 'assets', 'stl', 'threebellow_holder.stl');    
    stlPath2 =fullfile(currentDir, '..', 'assets', 'stl', 'threebellow_connector_res.stl');    

    % try
        obj1 = Gmodel(stlPath1,'Texture',matcap_egg,'Shading','Face');
        obj1 = obj1.render;

        obj2 = Gmodel(stlPath2,'Texture',matcap_egg,'Shading','Face');
        obj2 = obj2.render;

        if nargout > 2
            for ii = 1:(nargout-2)
                obj = Gmodel(file2Path,'Texture',matcap_egg,'Shading','Face');
                obj = obj.render;
                varargout(ii) = {obj};
                % 1
            end
        end
    % end
end
