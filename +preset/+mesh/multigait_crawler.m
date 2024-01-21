function msh = multigait_crawler(varargin)

    p = inputParser;
    addOptional(p,'n',1.5);
    parse(p,varargin{:});
    
    currentDir = fileparts(mfilename('fullpath'));
    imagePath = fullfile(currentDir, '..', 'assets', 'img', 'multigait_crawler.png');

    % try
        msh = Mesh(imagePath,'BdBox',[0,150,0,12],...
            'SimplificationTolerance',0.05,'ElementSize',p.Results.n,...
            'ElementOrder','linear');
        msh = msh.generate(); 
    % catch me
        
    % end
end