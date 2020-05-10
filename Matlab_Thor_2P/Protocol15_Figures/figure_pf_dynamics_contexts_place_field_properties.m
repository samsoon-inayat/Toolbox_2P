function figure_pf_dynamics_contexts_place_field_properties(fn,allRs,ccs)
adata = evalin('base','data');
mData = evalin('base','mData');
colors = mData.colors;
sigColor = mData.sigColor;
mData.axes_font_size = mData.axes_font_size-1;
selAnimals = 1:11;
mData.belt_length = adata{selAnimals(2)}{1}{1}.belt_length;
n = 0;

%%
trials = 3:10;
trials10 = 3:9;
for ii = 1:3%length(data)
    distDi = []; pwsi = []; pcsi = []; rsi = []; sii = [];
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
    end
    distD{ii} = distDi; pw{ii} = pwsi; pc{ii} = pcsi; Rsq{ii} = rsi; SI{ii} = sii;
end
pcs5 = distD; distD = [];
selContexts = [2 3];
selCells = selectCells15(selAnimals,'Common23');
% selCells = selectCells10(selAnimals,'Disrupted_C1');
% selCells = selectCells10(selAnimals,'Remained_C1');

%% Percentage of place cells vs place field centers

runThis = 0;
if runThis
distD = [];
for ii = 1:length(pc)
    tempV = pc{ii}(selCells);
    distD{ii} = tempV(~isnan(tempV));
end    
hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 7.5 6 2],'color','w');hold on;
[ha,hb,hca,sigR] = plotDistributions(distD,'colors',colors,'maxY',30,'cumPos',[0.6 0.4 0.23 0.35],'min',0,'incr',10,'max',mData.belt_length,'BaseValue',0.15);
set(gca,'TickDir','out','FontSize',mData.axes_font_size,'FontWeight','Bold');
hx = xlabel('Place Field Center (cm)');changePosition(hx,[0 0.0 0]);
hy = ylabel('Percentage (%)');changePosition(hy,[1 -1 0]);
legs = {'Context 1','Context 2','Context 3','Context 4',[60 5 30 3]};
% sigR = significanceTesting(distD);
putLegend(gca,legs,'colors',colors,'sigR',{sigR,'ks',sigColor,10});
changePosition(gca,[-0.04 0.09 0.11 -0.06]);
axes(hca);set(gca,'FontSize',10);
save_pdf(hf,mData.pdf_folder,'figure_percent_place_cells_vs_centers_contexts_15.pdf',600);
return;
end

%% scatter place field widths vs centers
runThis = 0;
if runThis
% hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 4 6 2],'color','w');hold on;
hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 4 3.5 2],'color','w');hold on;
ya = []; xa = []; gr = [];
for iii = 1:length(selContexts)
    ii = selContexts(iii);
    pwidths = pw{ii}(selCells)'; 
    pcenters = pc{ii}(selCells)';
%     nanspw = isnan(pwidths); nanspc = isnan(pcenters); nanspcs = (pcenters<3); nans = nanspw | nanspc | nanspcs; 
%     pwidths = pwidths(~nans); pcenters = pcenters(~nans);
    scatter(pcenters,pwidths,20,'.','MarkerEdgeColor',colors{ii},'MarkerFaceColor','none');
    ft = fittype('poly1');
    [ftc,gof,output] = fit(pcenters,pwidths,ft);    rsq(ii) = gof.rsquare;
    co = coeffvalues(ftc);
    pwf = ft(co(1),co(2),pcenters);
    hold on;plot(pcenters,pwf,'color',colors{ii},'LineWidth',1)
    xa = [xa;pcenters]; ya = [ya;pwidths]; gr = [gr;ones(size(pcenters))*iii];
end
sigR = do_ancova(xa,ya,gr);
ylim([0 95]);%xlim([0 mData.belt_length]);
set(gca,'TickDir','out','FontSize',mData.axes_font_size,'FontWeight','Bold');
hx = xlabel('Place Field Center (cm)');changePosition(hx,[0 3.0 0]);
hy = ylabel('Place Field Width (cm)');changePosition(hy,[1 -14.0 0]);
for iii = 1:length(selContexts)
    ii = selContexts(iii);
    legs{ii} = sprintf('C %d',ii);%{'Context 1','Context 2','Context 3','Context 4',[30 3 35 3]};
end
legs = legs(selContexts);
legs{length(legs)+1} = [35 5 80 7];putLegend(gca,legs,'colors',colors(selContexts),'sigR',{sigR,'ancovas',sigColor,9}); text(30,90,'Slope','FontWeight','Bold');
changePosition(gca,[0.02 0.09 0.01 -0.06]);%changePosition(gca,[-0.03 0.09 0.09 -0.06]);
save_pdf(hf,mData.pdf_folder,'figure_scatter_centers_widths_contexts_15.pdf',600);
return;
end

%% scatter zMI vs centers
runThis = 1;
if runThis
% hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 4 6 2],'color','w');hold on;
hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 4 3.5 2],'color','w');hold on;
ya = []; xa = []; gr = [];
for iii = 1:length(selContexts)
    ii = selContexts(iii);
    SIs = SI{ii}(selCells)';
    pcenters = pc{ii}(selCells)';
    scatter(pcenters,SIs,20,'.','MarkerEdgeColor',colors{ii},'MarkerFaceColor','none');
    ft = fittype('poly1');
    [ftc,gof,output] = fit(pcenters,SIs,ft);  rsq(ii) = gof.rsquare;
    co = coeffvalues(ftc);
    pwf = ft(co(1),co(2),pcenters);
    hold on;plot(pcenters,pwf,'color',colors{ii},'LineWidth',1)
    xa = [xa;pcenters]; ya = [ya;SIs]; gr = [gr;ones(size(pcenters))*iii];
end
sigR = do_ancova(xa,ya,gr);
ylim([0 45]); %xlim([0 100]);
set(gca,'TickDir','out','FontSize',mData.axes_font_size,'FontWeight','Bold');
hx = xlabel('Place Field Center (cm)');changePosition(hx,[0 1.0 0]);
hy = ylabel('Mutual Info. (Z Score)');changePosition(hy,[0 -5.0 0]);
for iii = 1:length(selContexts)
    ii = selContexts(iii);
    legs{ii} = sprintf('C %d',ii);%{'Context 1','Context 2','Context 3','Context 4',[30 3 35 3]};
end
legs = legs(selContexts);
legs{ii+1} = [70 5 40 4];putLegend(gca,legs,'colors',colors(selContexts),'sigR',{sigR,'ancovai',sigColor,9}); text(65,45,'Intercept','FontWeight','Bold'); 
% legs{end} = [70 5 33 2];putLegend(gca,legs,'colors',colors,'sigR',{sigR,'ancovas',sigColor,9}); text(65,36,'Slope','FontWeight','Bold');

changePosition(gca,[0.02 0.09 0.03 -0.06]);%changePosition(gca,[-0.05 0.09 0.11 -0.06]);
save_pdf(hf,mData.pdf_folder,'figure_scatter_centers_zMI_contexts_15.pdf',600);
return;
end


