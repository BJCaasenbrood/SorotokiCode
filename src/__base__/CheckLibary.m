function LatestRelease = CheckLibary(lib,ver)
ver_ = vernum(lib);
% if length(lib) > 7
% cprintf('Keyword',['* \t -',lib,spc,'version: ',ver,spc]);
% else
cout('Keyword',['* \t -',lib,spc,'version: ',ver,spc]);
% end

if strcmp(ver,ver_)
cout('Green', '[ok]') 
LatestRelease = true;
else
cprintf('UnterminatedStrings', 'outdated...')
LatestRelease = false;
end
cprintf('Text','\n')
end

%-------------------------------------------------------------------- SPACE
function x = spc
x = '\t | ';
end

