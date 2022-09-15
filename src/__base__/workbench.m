function workbench()
list = {'Go to workbench',...
        'Add to workbench',...
        'Remove workbench'};
    
% pth = userpath;
% fopen([pth,'workplaces.m']);
jobs = workplaces;    
h = action(jobs);
cd(jobs{h});
end

