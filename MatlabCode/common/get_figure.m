function hf = get_figure(fn,pos)

hf = figure(fn);clf;
set(gcf,'Units','Inches');
set(gcf,'Position',pos,'color','w');
hold on;