function msh = suzumori_walker(varargin)

    p = inputParser;
    addOptional(p,'n',1500);
    parse(p,varargin{:});
    
    currentDir = fileparts(mfilename('fullpath'));
    imagePath = fullfile(currentDir, '..', 'assets', 'suzumori_walker.png');

    try
        msh = Mesh(imagePath,'BdBox',[0,80,0,17],...
            'SimplificationTolerance',0.025,'NElem',p.Results.n,...
            'ElementOrder','Quadratic');
        msh = msh.generate(); 
    end
end