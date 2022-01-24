function folderFound = CheckLibary(lib,fnc)
folderFound = true;
cout('Keyword',['* \t -',lib,spc,'version: ',soropatch(1),spc]);

try 
    fnc(1);
catch
    folderFound = false;
end

if folderFound
    cout('Green', '[ok]') 
else
    cout('UnterminatedStrings', 'outdated...')
end

cout('Text','\n')
end

%-------------------------------------------------------------------- SPACE
function x = spc
x = '\t | ';
end

