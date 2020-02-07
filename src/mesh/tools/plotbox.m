function h = plotbox
h = gca;
x = h.XLim;
y = h.YLim;

hold on;
h = plot([x(1),x(2),x(2),x(1),x(1)],[y(1),y(1),y(2),y(2),y(1)],'k-',...
    'linewidth',1.5);

end

