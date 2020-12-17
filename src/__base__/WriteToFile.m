function WriteToFile(str)
global FID;
if ismac
    %
elseif isunix
    str = strrep(str,'\','/');
elseif ispc
    str = strrep(str,'/','\');
else
    %
end
fprintf(FID, ['addpath(genpath(''',str, '''))\n']);
pause(0.01);
cprintf('Keyword',[str,'\n']);
end
