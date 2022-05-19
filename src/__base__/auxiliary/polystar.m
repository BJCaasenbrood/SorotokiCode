function P = polystar(M,r1,r2)

if nargin < 2
   r1 = 0.333*sqrt(2);
   r2 = 1;
end

if nargin < 3
   r2 = r1;
   r1 = 0.333*sqrt(2)*r2;
end

N = 2*pi/(M);
angles1 = 0:N:2*pi;
angles2 = (N/2:N:2*pi+N/2);
x1 = r1 * cos(angles1 - .5*pi);
y1 = r1 * sin(angles1 - .5*pi);
x2 = r2 * cos(angles2 - .5*pi);
y2 = r2 * sin(angles2 - .5*pi);

P = [(reshape([x1;x2], 1, []))',...
     (reshape([y1;y2], 1, []))'];
 
P = flipud(P); 
P(end,:) = [];
end

