function [x,y] = recoverXY()
fig = gcf;
axObjs = fig.Children;
dataObjs = axObjs.Children;
x = dataObjs(1).XData;
y = dataObjs(1).YData;
end