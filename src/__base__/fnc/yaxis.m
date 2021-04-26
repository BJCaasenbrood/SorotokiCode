function yaxis(Label,varargin)
if isempty(varargin)
    ylabel(Label,'interpreter','latex','fontsize',14);
else
    ylabel([Label, ' ','(',varargin{1},')'],'interpreter','latex','fontsize',14);
end
end