function output = action(list)

CommandSymbol = '*';

fprintf('* Please, make a selection: \n');

for i = 1:length(list)
   cprintf('Green',['\t','  ','[',num2str(i),'] ',list{i}]);
   fprintf('\n'); 
end

fprintf(CommandSymbol);
fprintf('  '); 

cprintf('Text','>> ');
output = input('');
end


