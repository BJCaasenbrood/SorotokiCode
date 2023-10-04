function out = runSorotokiTest(prompt)
flag = [];
base = 'src/Sorotoki';
testsuite = {'Sdf'; 'Mesh'; 'Fem'};

% eliminate suites not in prompt
if ~isempty( prompt )
    redux = ismember(lower(testsuite), prompt);
    testsuite = testsuite(redux);
end

% loop over testsuites
for ii = 1:numel(testsuite)
    flag = navigateAndTest(base, testsuite{ii}, flag);
end

if all(flag == 1)
    fprintf('All tests passed!\n ');
    out = true;
else
    out = false;
end

end

%% nagivates to library and tests
function flag = navigateAndTest(base, lib, flag)
    sorotoki cd; 
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
