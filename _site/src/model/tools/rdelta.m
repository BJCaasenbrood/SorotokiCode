function y = rdelta(x,h)
y = zeros(numel(x),1);
for ii = 1:numel(x)
   if abs(x(ii)) < h/2, y(ii) = 1/h;
   else, y(ii) = 0;
   end
end

end

