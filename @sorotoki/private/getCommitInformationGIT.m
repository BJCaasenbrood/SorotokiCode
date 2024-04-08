function [sha, date] = getCommitInformationGIT(repo)

    api   = 'https://api.github.com/repos/bjcaasenbrood/';
    where = '/git/refs/heads/master';

    % assert(ismember(repo,soroPackages), 'Repo not a member from Sorotoki Packages');
    responseMaster = webread([api, repo, where]);
    response = webread(responseMaster.object.url);

    reply = response.author;
    sha = response.sha;
    date = datetime(reply.date, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');

end