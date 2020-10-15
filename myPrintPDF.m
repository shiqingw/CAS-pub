function []=myPrintPDF(fig_h,leg_h,filename)
% Prints a PDF of a figure, with the correct Fonts/LineWidth etc.
% h = handle of the figure


 
hAxes=get(fig_h,'CurrentAxes');
set(get(hAxes,'xlabel'),'Interpreter','Latex','Fontsize',14)
set(get(hAxes,'ylabel'),'Interpreter','Latex','Fontsize',14)
set(get(hAxes,'title'),'Interpreter','Latex','Fontsize',14)
if ~isempty(leg_h)
    set(leg_h,'Interpreter','Latex','Fontsize',14)
end
% set(gco,'Interpreter','Latex','Fontsize',14)
%   set(hAxes,...
% 'Units','normalized',...
% 'Position',[.1 .1 .85 .85],'xLimMode','auto','yLimMode',...
% 'auto','zLimMode','auto','XGrid','on','YGrid','on','ZGrid','on')
plots=get(hAxes,'Children');
set(plots,'LineWidth',3)
set(fig_h,'Units','inches',...
    'Position',[1 3 7 5],...
    'PaperPositionMode','auto','PaperUnits','inches','PaperSize',[7 5])
print('-dpdf',filename,'-r90')