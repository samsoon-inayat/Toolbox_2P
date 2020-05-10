thisCols = {'r','b'};
xs = [1 2 3 4 10 11];
a_171034 = [25.8 25.8 25.9 25.9 26.7 26.9];
a_170865 = [31.4 31.2 31.1 31.1 31.5 31.0];

pp.numberOfRows = 1;
pp.spaceBetweenRows = [0.05];
pp.rowHeights = [0.6];

pp.numberOfCols = [1];
pp.colWidths = {[0.75]};
pp.spaceBetweenCols = {[0]};

pp.leftOffset = 0.17;
pp.bottomOffset = 0.25;

pp = getPanelsProps(pp);


hf = makeFigureWindow(1,[5 5 3.5 2],1);

pp = makeAxes(hf,pp,1,1,[0 0 0 0]);

plot(xs,a_171034,'marker','o','linewidth',1.5,'color',thisCols{1});hold on;
plot(xs,a_170865,'marker','o','linewidth',1.5,'color',thisCols{2});
xlabel('Day');
ylabel('Weight (g)');
ylim([22 37]);
xlim([0 12])

x1 = 5; x2 = x1+0.75; y1 = 34:1.5:37; y2 = y1;
legs = {'Animal 2','Animal 1'};
for ii = 2:-1:1
    plot([x1 x2],[y1(ii) y2(ii)],'color',thisCols{ii});
    text(x2+0.2,y1(ii),legs{ii},'Color',thisCols{ii},'FontSize',9);
end
box off;
set(gca,'TickDir','out','FontSize',10,'FontWeight','Bold')
title('Weight vs Days')
saveFigAsPDF(hf,'weightVsDays.pdf',0);