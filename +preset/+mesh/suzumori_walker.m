function msh = suzumori_walker(varargin)

    p = inputParser;
    addOptional(p,'n',1500);
    parse(p,varargin{:});
    
    currentDir = fileparts(mfilename('fullpath'));
    imagePath = fullfile(currentDir, '..', 'assets', 'img','suzumori_walker.png');

    msh = Mesh(imagePath,'BdBox',[0,80,0,17],...
        'SimplificationTolerance',0.025,'NElem',p.Results.n,...
        'ElementOrder','linear');
    
    msh = msh.generate();
end