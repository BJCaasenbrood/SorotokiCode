function checkToolboxes(reqToolboxes)

    v = ver;
    installedToolboxes = {v.Name};
    missingToolboxes   = {};
    
    for i = 1:numel(reqToolboxes)
        if ~ismember(reqToolboxes{i},installedToolboxes)
            if isempty(missingToolboxes)
                disp(msg('missing_toolboxes'));
            end
            disp([' - ', reqToolboxes{i}]);
            missingToolboxes{end+1} = reqToolboxes{i};
        end
    end
    
    if ~isempty(missingToolboxes)
        assert(false, msg('request_to_install_toolbox'))
    end

end