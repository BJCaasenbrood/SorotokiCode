function print(str)
if isa(str,'string'), cout(strcat(">> ",str));
elseif isa(str,'char'), cout(['>> ',str]);
end
end

