function msh = sorotoki_logo(varargin)
% PNEUNET Generates a mesh using the Mesh class.

    p = inputParser;

    p.addOptional('n',1000); % number of elements
    parse(p,varargin{:});
    
    currentDir = fileparts(mfilename('fullpath'));
    imagePath = fullfile(currentDir, '..', 'assets', ...
        'img','sorotoki_bwlogo.png');

    try
        msh = Mesh(imagePath,'BdBox',[0,120,0,120],...
            'NElem',p.Results.n,'ElementOrder','linear');
        msh = msh.generate();
    end
end
