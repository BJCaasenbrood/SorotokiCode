function [dataObjs] = recoverXY(varargin)
fig = gcf;
axObjs = fig.Children;
dataObjs = axObjs.Children;
% x = dataObjs(1).XData;
% y = dataObjs(1).YData;

if ~isempty(varargin)
    if isa(varargin{1},'char')
        csvwrite(varargin{1},...
            [dataObjs(1).XData(:) dataObjs(1).YData(:)]);
    end
end

end