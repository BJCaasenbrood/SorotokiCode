function LatestRelease = CheckLibary(lib,ver)
ver_ = vernum(lib);
if length(lib) > 7
cprintf('Keyword',['* \t -',lib,spc,'version: ',ver,spc]);
else
cprintf('Keyword',['* \t -',lib,'\t',spc,'version: ',ver,spc]);
end

if strcmp(ver,ver_)
cprintf('Green', '[ok]') 
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

