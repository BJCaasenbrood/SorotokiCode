function workbench()
list = {'Go to workbench',...
        'Add to workbench',...
        'Remove workbench'};
    
jobs = workplaces;    
h = action(jobs);
cd(jobs{h});
end

