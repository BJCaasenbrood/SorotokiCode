function id = whoClasses(ClassName)
evalin('base', 'clearvars ans;');
s = evalin('base', 'whos');
data_collector = zeros(size(s));
id = cell(1);
j = 1;
for i = 1:length(s)
    k = s(i).class;
    data_collector(i) = strcmp(k,ClassName); 
    
    if data_collector(i) == 1
        id{j} = evalin('base', (s(i).name));
        j = j+1;
    end
end

end

