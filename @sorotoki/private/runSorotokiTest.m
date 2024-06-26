function out = runSorotokiTest(prompt, installPath, base, testsuite)
flag = [];

% eliminate suites not in prompt
try 
    if ~isempty( prompt{1} )
        redux = ismember(lower(testsuite), prompt);
        testsuite = testsuite(redux);
    end
catch
    % ...
end

% loop over testsuites
warning off;
for ii = 1:numel(testsuite)
    cd(installPath);
    flag = navigateAndTest(base, testsuite{ii}, flag);
end
warning on;

if all(flag == 1)
    fprintf('All tests passed!\n');
    out = true;
else
    out = false;
end

end

%% nagivates to library and tests
function flag = navigateAndTest(base, lib, flag)
    cd([base,lib]);
    fprintf(['Running ', lib, 'Test ']);
    test = runtests(pwd, 'OutputDetail',1);
    flag(end+1) = showTestResults(test);    
end

%% processes test results
function isPassed = showTestResults(test)
    log = Log;
    log.hline();
    name = {test.Name};

    flag = ones(numel(name),1); 
    flag(vertcat(test.Failed)) = 0;

    for i = 1:numel(name)
        log.bool(name{i},flag(i),{'Passed','Failed'});
        pause(1e-2);
    end
    log.hline();

    % assertSuccess(test);
    isPassed = logical(cumprod(flag));
    isPassed = isPassed(end);
end