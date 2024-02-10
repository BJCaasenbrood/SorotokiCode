function msh = pneunet_AMOLF(varargin)
% PNEUNET Generates a mesh using the Mesh class.

    p = inputParser;

    p.addOptional('n',2); % number of elements
    parse(p,varargin{:});
    
    currentDir = fileparts(mfilename('fullpath'));
    imagePath = fullfile(currentDir, '..', 'assets', 'img', 'Amolf_Actuator_25x118.png');

    try
        msh = Mesh(imagePath,'BdBox',[0,118,0,25],...
            'ElementSize',p.Results.n,'ElementOrder','linear','SimplificationTolerance',0.1);
        msh = msh.generate();
    end
end
