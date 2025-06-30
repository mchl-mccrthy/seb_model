% Format figures
%
% Michael McCarthy 2025

function formatfigure(width,height)

box on; set(gca,'Layer','top') % Plot box on top
set(gca,'linewidth',0.5); set(gca,'fontsize',10)
set(gcf,'units','centimeters','position',[0 0 width height]); hold off

end