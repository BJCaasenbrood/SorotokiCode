clr;

C = 3;
M = 4;
Y = affinechebyspace(1e3,M,C,false);

% plotting
for ii = 1:size(Y,2)
   plot(Y(:,ii),'LineW',1.5); hold on;
end

sorocolor;