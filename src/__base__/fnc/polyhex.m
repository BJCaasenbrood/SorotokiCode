function P = polyhex(M,r1)

if nargin < 2
   r1 = 1;
end

N = 2*pi/(M);
angles1 = 0:N:2*pi;
x1 = r1 * cos(angles1 - .5*pi);
y1 = r1 * sin(angles1 - .5*pi);

P = [(reshape([x1], 1, []))',...
     (reshape([y1], 1, []))'];
 
P = flipud(P); 
%P(end,:) = [];
end

