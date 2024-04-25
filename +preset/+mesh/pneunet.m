function msh = pneunet(varargin)
% PNEUNET Generates a mesh using the Mesh class.

    p = inputParser;
    p.addOptional('n',2); % number of elements
    parse(p,varargin{:});
    
    currentDir = fileparts(mfilename('fullpath'));
    imagePath = fullfile(currentDir, '..', 'assets', 'img', 'pneunet.png');

    try
        msh = Mesh(imagePath,'BdBox',[0,120,0,20],...
            'ElementSize',p.Results.n,'ElementOrder','linear','SimplificationTolerance',0.0125);
        msh = msh.generate();
    end
end
