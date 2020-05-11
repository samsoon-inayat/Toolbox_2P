function figure_place_field_properties(fn,allRs,ccs)
adata = evalin('base','data');
mData = evalin('base','mData');
colors = mData.colors;
sigColor = mData.sigColor;
% selAnimals = 1:3;
% selAnimals = 6:11;
selAnimals = 1:8;
mData.belt_length = adata{selAnimals(1)}{1}{1}.belt_length;
n = 0;

%%
trials = 3:10;
trials10 = 3:9;
for ii = 1:4%length(data)
    distDi = []; pwsi = []; pcsi = []; rsi = []; sii = []; hafdii = []; hifdii = []; 
    for jj = 1:length(selAnimals)
        [ii jj selAnimals(jj)]
        [tempD cnsjj] = getVariableValues(adata{selAnimals(jj)},'placeCells5',ii);
        distDi = [distDi tempD];
        [tempD cnsjj] = getVariableValues(adata{selAnimals(jj)},'pws',ii);
        pwsi = [pwsi tempD];
        [tempD cnsjj] = getVariableValues(adata{selAnimals(jj)},'centers',ii);
        pcsi = [pcsi tempD];
        [tempD cnsjj] = getVariableValues(adata{selAnimals(jj)},'rs',ii);
        rsi = [rsi tempD];
        [tempD cnsjj] = getVariableValues(adata{selAnimals(jj)},'SI',ii);
        sii = [sii tempD];
        [tempD cnsjj] = getVariableValues(adata{selAnimals(jj)},'HaFD',ii);
        hafdii = [hafdii tempD];
        [tempD cnsjj] = getVariableValues(adata{selAnimals(jj)},'HiFD',ii);
        hifdii = [hifdii tempD];
    end
    distD{ii} = distDi; pw{ii} = pwsi; pc{ii} = pcsi; Rsq{ii} = rsi; SI{ii} = sii; HaFD{ii} = hafdii;  HiFD{ii} = hifdii;
end
pcs5 = distD; distD = [];

%% Percentage of diff zMI
runThis = 1;
if runThis
distD = [];
combs = nchoosek(1:4,2);
for ii = 1:size(combs,1)
    first = combs(ii,1); second = combs(ii,2);
    tempV = SI{second} - SI{first};
    distD{ii} = tempV(~isnan(tempV));
end    
hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 7.5 6 2.5],'color','w');hold on;
[ha,hb,hca,sigR] = plotDistributions(distD,'colors',colors,'maxY',50,'cumPos',[0.6 0.4 0.23 0.35],'min',-10,'incr',2,'max',20,'BaseValue',0.15);
set(gca,'TickDir','out','FontSize',mData.axes_font_size,'FontWeight','Bold');
hx = xlabel('Difference zMI');changePosition(hx,[0 0.0 0]);
hy = ylabel('Percentage (%)');changePosition(hy,[0.25 -1 0]);
combs
legs = {'C12','C13','C14','C23','C24','C34',[6 1 50 3]};
% sigR = significanceTesting(distD);
putLegend(gca,legs,'colors',colors,'sigR',{sigR,'anova',sigColor,10});
changePosition(gca,[-0.04 0.09 0.11 -0.06]);
axes(hca);set(gca,'FontSize',10);
save_pdf(hf,mData.pdf_folder,'figure_percent_cells_vs_diff_zMI.pdf',600);
return;
end

%% Percentage of diff Rs
runThis = 1;
if runThis
distD = [];
combs = nchoosek(1:4,2);
for ii = 1:size(combs,1)
    first = combs(ii,1); second = combs(ii,2);
    tempV = Rsq{second} - Rsq{first};
    distD{ii} = tempV(~isnan(tempV));
end    
hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 7.5 6 2.5],'color','w');hold on;
[ha,hb,hca,sigR] = plotDistributions(distD,'colors',colors,'maxY',50,'cumPos',[0.6 0.4 0.23 0.35],'min',-1.5,'incr',0.2,'max',1.5,'BaseValue',0.15);
set(gca,'TickDir','out','FontSize',mData.axes_font_size,'FontWeight','Bold');
hx = xlabel('Difference Rsq');changePosition(hx,[0 0.0 0]);
hy = ylabel('Percentage (%)');changePosition(hy,[0 -1 0]);
combs
legs = {'C12','C13','C14','C23','C24','C34',[0.1 0.1 48 3]};
% sigR = significanceTesting(distD);
putLegend(gca,legs,'colors',colors,'sigR',{sigR,'ks',sigColor,10});
changePosition(gca,[-0.04 0.09 0.11 -0.06]);
axes(hca);set(gca,'FontSize',10);
save_pdf(hf,mData.pdf_folder,'figure_percent_cells_vs_diff_zMI.pdf',600);
return;
end

