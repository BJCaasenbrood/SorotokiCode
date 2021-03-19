function y = zci(x)
% ZCI finds zero crossing index in a vector
y = find(x(:).*circshift(x(:), [-1 0]) <= 0);
end

