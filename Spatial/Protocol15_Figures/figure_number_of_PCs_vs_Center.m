function figure_number_of_PCs_vs_Center(fn,allRs,ccs)
adata = evalin('base','data');
mData = evalin('base','mData');
colors = mData.colors;
sigColor = mData.sigColor;
% selAnimals = 1:11;
selAnimals = [2 4 6 8 10];
mData.belt_length = adata{selAnimals(2)}{1}{1}.belt_length;
n = 0;

%% collect values
for jj = 1:length(selAnimals)
    for ii = 1:3%length(data)
        [ii jj selAnimals(jj)]
        [distDi cnsjj] = getVariableValues(adata{selAnimals(jj)},'placeCells5',ii);
        [pwsi cnsjj] = getVariableValues(adata{selAnimals(jj)},'pws',ii);
        [pcsi cnsjj] = getVariableValues(adata{selAnimals(jj)},'centers',ii);
        [rsi cnsjj] = getVariableValues(adata{selAnimals(jj)},'rs',ii);
        [sii cnsjj] = getVariableValues(adata{selAnimals(jj)},'SI',ii);
        pc5{jj,ii} = distDi; pw{jj,ii} = pwsi; pc{jj,ii} = pcsi; Rsq{jj,ii} = rsi; SI{jj,ii} = sii;
    end
end

for jj = 1:length(selAnimals)
    allpcsU = pc5{jj,1}; 
    for ii = 2:3%length(data)
         allpcsU = allpcsU | pc5{jj,ii}; 
         lastPCs = pc5{jj,ii-1}; currentPCs = pc5{jj,ii};
    end
    percent_pcs(jj,1) = 100*sum(allpcsU)/length(pc5{jj,1});
%     percent_pcsi(jj,:) = 100*pcs(jj,:)./length(pc5{jj,1});
    upcs(jj,1) = sum(allpcsU);
end

for jj = 1:length(selAnimals)
    for ii = 2:3%length(data)
         lastPCs = pc5{jj,ii-1}; currentPCs = pc5{jj,ii};
         remained(jj,ii) = 100*sum(lastPCs & currentPCs)/sum(currentPCs);
         disrupted(jj,ii) = 100*sum(lastPCs & ~currentPCs)/sum(lastPCs);
         newones(jj,ii) = 100*sum(~lastPCs & currentPCs)/sum(currentPCs);
         newPCs{jj,ii-1} = ~lastPCs & currentPCs;
    end
end

%%
runThis = 0;
if runThis
    space_bins = 20:20:170;
    for jj = 1:length(selAnimals)
        for ii = 1:2
            barsd{ii}(jj,:) = hist(pc{jj,ii}(find(newPCs{jj,ii})),space_bins);
        end
%         bars2(jj,:) = hist(pc{jj,2},space_bins);
%         bars3(jj,:) = hist(pc{jj,3},space_bins);
    end
    hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 5 3.5 2.5],'color','w');
    hold on;
    for ii = 1:2
        d = mean(barsd{ii});
        s = std(barsd{ii})/sqrt(length(selAnimals));
        shadedErrorBar(space_bins,d,s,colors{ii},0.5);
    end
    xlim([10 150]);
    return;
end
n = 0;

%% Percentage of place cells vs place field centers
runThis = 1;
if runThis
distD = [];
for ii = 1:2
    tempV = pc{ii}(find(newPCs{ii}));
    distD{ii} = tempV(~isnan(tempV));
end    
hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 7.5 6 2.5],'color','w');hold on;
[ha,hb,hca,sigR] = plotDistributions(distD,'colors',colors,'maxY',30,'cumPos',[0.6 0.4 0.23 0.35],'min',0,'incr',10,'max',mData.belt_length,'BaseValue',0.15);
set(gca,'TickDir','out','FontSize',mData.axes_font_size,'FontWeight','Bold');
hx = xlabel('Place Field Center (cm)');changePosition(hx,[0 0.0 0]);
hy = ylabel('Percentage (%)');changePosition(hy,[1 -1 0]);
legs = {'Context 1','Context 2','Context 3',[60 5 30 3]};
% sigR = significanceTesting(distD);
putLegend(gca,legs,'colors',colors,'sigR',{sigR,'ks',sigColor,10});
changePosition(gca,[-0.04 0.09 0.11 -0.06]);
axes(hca);set(gca,'FontSize',10);
save_pdf(hf,mData.pdf_folder,'figure_percent_place_cells_vs_centers_15.pdf',600);
return;
end

%%
runThis = 1;
if runThis
    
    return;
end
n = 0;
