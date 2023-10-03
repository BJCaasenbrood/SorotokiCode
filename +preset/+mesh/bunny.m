function msh = bunny(varargin)

    p = inputParser;
    addOptional(p,'n',1000);
    parse(p,varargin{:});
    
    currentDir = fileparts(mfilename('fullpath'));
    modelPath = fullfile(currentDir, '..', 'assets', 'stl', 'Bunny.stl');

    try
        msh = Mesh(modelPath,'Hmesh',[5,10]);
        msh = msh.generate;
    catch me
        
    end
end