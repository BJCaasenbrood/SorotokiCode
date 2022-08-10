function f = fig(varargin)
n = 101;
B = [10,9.58];

for ii = 1:numel(varargin)
   if numel(varargin{ii}) == 1
       n = varargin{ii};
   elseif numel(varargin{ii}) == 2
       B = round(varargin{ii},2);
       if B(1) > 14
           B = log(B)/log(2);
       end
   end
end

f = figure(n);

f.Units = 'pixel';
Dx = 2^B(1);
Dy = 2^B(2);
b = 0.5715;
a = 512;
f.Position = round([a,a,b*Dx,b*Dy]);
set(gca,'LineW',1.5);
background();
hold on;

end

