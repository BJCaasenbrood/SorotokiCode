function obj = tennis_racket

    currentDir = fileparts(mfilename('fullpath'));
    % stlPath = fullfile(currentDir, 'tennis-racket-theorem.stl');
    stlPath =fullfile(currentDir, '..', 'assets', 'stl', 'tennis-racket-theorem.stl');

    obj = Gmodel(stlPath,'Texture', matcap_diffuse(0.5),'Shading','Face');
    obj = obj.render;
end