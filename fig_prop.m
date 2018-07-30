function fig_prop(width,height)

figure('units','normalized','outerposition',[0 0 1 1]);
fig = gcf;
set(fig, 'Units', 'inches')
fig.PaperUnits = 'inches';
set(fig,'PaperPositionMode','manual');
fig.PaperPosition = [0 0 width height];
fig.PaperSize = [width height];

set(gcf,'renderer','painters');

end