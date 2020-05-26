function output = action(list)

CommandSymbol = '*';

CallDisplay('please, make a selection:');

for i = 1:length(list)
   cprintf('Strings',['\t','  ','[',num2str(i),'] ',list{i}]);
   fprintf('\n'); 
end

fprintf(CommandSymbol);
fprintf('*    '); 

cprintf('Text','input: ');
output = input('');
end


