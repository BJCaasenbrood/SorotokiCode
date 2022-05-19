function xaxis(Label,varargin)
if isempty(varargin)
    xlabel(Label,'interpreter','latex','fontsize',14);
else
    xlabel([Label, ' ','(',varargin{1},')'],'interpreter','latex','fontsize',14);
end
end