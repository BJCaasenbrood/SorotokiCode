function msh = diamond_bot(varargin)

    p = inputParser;
    addOptional(p,'n',3.5);
    parse(p,varargin{:});
    
    currentDir = fileparts(mfilename('fullpath'));
    imagePath = fullfile(currentDir, '..', 'assets', 'img', 'diamond_bot.png');

    msh = Mesh(imagePath,'BdBox',[0,170,0,100],...
        'SimplificationTolerance',0.025,'ElementSize',p.Results.n,...
        'ElementOrder','linear');
    msh = msh.generate(); 
end