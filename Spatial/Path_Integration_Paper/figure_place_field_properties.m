function figure_place_field_properties(fn,allRs,ccs)
ei = evalin('base','ei10');
mData = evalin('base','mData');
colors = mData.colors;
sigColor = mData.sigColor;
selAnimals = [1:4 9];
n = 0;
%%
runthis = 1;
if runthis
%     cN = 1;
    planeNumbers = 'All';
    maxDistTime = [142 15];
    contextNumber = ones(1,4).*1;
    
    stimMarkers = {'air'};%,'air','belt','airI'};
    stimLabel = {'Air','Air','Belt','InterTrial'};
    xstimLabel = {'AirD','AirT','BeltD','AirIT'};
    rasterTypes = {'dist','time','dist','time'};
    rasterLabel = {'Distance','Time','Distance','Time'};
    selCells = 'areCells';
    allVals = []; gAllVals = [];
    all_data = [];
    for cnsi = 1:4
        for ss = 1:length(stimMarkers)
            distD = [];
            for jj = 1:length(selAnimals)
                [pcs cns areCells] = getParamValues('placeCells3',ei(selAnimals(jj)),planeNumbers,cnsi,'belt','dist',selCells,maxDistTime);
%                 [clus] = getParamValues('cluster3',ei(selAnimals(jj)),planeNumbers,cnsi,'air','dist',selCells,maxDistTime);
                npcs(jj,cnsi) = sum(pcs)/length(pcs);
                [pws cns ACs] = getParamValues('place_field_properties.pws',ei(selAnimals(jj)),planeNumbers,cnsi,stimMarkers{ss},rasterTypes{ss},selCells,maxDistTime);
                [centers cns ACs] = getParamValues('place_field_properties.centers',ei(selAnimals(jj)),planeNumbers,cnsi,stimMarkers{ss},rasterTypes{ss},selCells,maxDistTime);
                [rs cns ACs] = getParamValues('place_field_properties.rs',ei(selAnimals(jj)),planeNumbers,cnsi,stimMarkers{ss},rasterTypes{ss},selCells,maxDistTime);
                centersI = centers > 3 & centers < 140;
                pwsI = pws < 120;
                inds = pcs & pwsI & centersI & rs>0.4;
                data{jj,cnsi} = [pws(inds) centers(inds)];
            end
        end
    end    
end

for ii = 1:size(data,2)
    all_pws = []; all_cs = [];
    for jj = 1:size(data,1)
        tpw = data{jj,ii}(:,1);
        all_pws = [all_pws;tpw];
        tpw = data{jj,ii}(:,2);
        all_cs = [all_cs;tpw]
    end
    all_pws_c{ii} = all_pws;
    all_cs_c{ii} = all_cs;
end
n=0;
%% scatter place field widths vs centers
runThis = 1;
if runThis
hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 4 6 2.5],'color','w');hold on;
ya = []; xa = []; gr = [];
for ii = 1:length(all_pws_c)
    pwidths = all_pws_c{ii};
    pcenters = all_cs_c{ii};
    scatter(pcenters,pwidths,20,'.','MarkerEdgeColor',colors{ii},'MarkerFaceColor','none');
    ft = fittype('poly1');
    [ftc,gof,output] = fit(pcenters,pwidths,ft);    rsq(ii) = gof.rsquare;
    co = coeffvalues(ftc);
    pwf = ft(co(1),co(2),pcenters);
    hold on;plot(pcenters,pwf,'color',colors{ii},'LineWidth',1)
    xa = [xa;pcenters]; ya = [ya;pwidths]; gr = [gr;ones(size(pcenters))*ii];
end
sigR = do_ancova(xa,ya,gr);
ylim([0 100]);
set(gca,'TickDir','out','FontSize',mData.axes_font_size,'FontWeight','Bold');
hx = xlabel('Place Field Center (cm)');changePosition(hx,[0 0.0 0]);
hy = ylabel('Place Field Width (cm)');changePosition(hy,[2.2 -3.0 0]);
for ii = 1:length(all_pws_c)
    legs{ii} = sprintf('C %d',ii);%{'Context 1','Context 2','Context 3','Context 4',[30 3 35 3]};
end
legs{ii+1} = [35 5 93 7];putLegend(gca,legs,'colors',colors,'sigR',{sigR,'ancovas',sigColor,9}); text(30,100,'Slope','FontWeight','Bold');
legs{ii+1} = [100 5 93 7];putLegend(gca,legs,'colors',colors,'sigR',{sigR,'ancovai',sigColor,9}); text(95,100,'Intercept','FontWeight','Bold');
changePosition(gca,[-0.03 0.09 0.09 -0.06]);
save_pdf(hf,mData.pdf_folder,'figure_scatter_centers_widths.pdf',600);
return;
end

%% scatter zMI vs centers
runThis = 0;
if runThis
hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[10 4 6 2.5],'color','w');hold on;
ya = []; xa = []; gr = [];
for ii = 1:length(pw)
    SIs = SI{ii}(find(pcs5{ii}))';
    pcenters = pc{ii}(find(pcs5{ii}))';
    scatter(pcenters,SIs,20,'.','MarkerEdgeColor',colors{ii},'MarkerFaceColor','none');
    ft = fittype('poly1');
    [ftc,gof,output] = fit(pcenters,SIs,ft);  rsq(ii) = gof.rsquare;
    co = coeffvalues(ftc);
    pwf = ft(co(1),co(2),pcenters);
    hold on;plot(pcenters,pwf,'color',colors{ii},'LineWidth',1)
    xa = [xa;pcenters]; ya = [ya;SIs]; gr = [gr;ones(size(pcenters))*ii];
end
sigR = do_ancova(xa,ya,gr);
ylim([0 35]);
set(gca,'TickDir','out','FontSize',mData.axes_font_size,'FontWeight','Bold');
hx = xlabel('Place Field Center (cm)');changePosition(hx,[0 0.0 0]);
hy = ylabel('Mutual Info. (Z Score)');changePosition(hy,[0 0.0 0]);
for ii = 1:length(pw)
    legs{ii} = sprintf('C %d',ii);%{'Context 1','Context 2','Context 3','Context 4',[30 3 35 3]};
end
legs{ii+1} = [130 5 33 2];putLegend(gca,legs,'colors',colors,'sigR',{sigR,'ancovai',sigColor,9}); text(125,36,'Intercept','FontWeight','Bold'); 
legs{end} = [70 5 33 2];putLegend(gca,legs,'colors',colors,'sigR',{sigR,'ancovas',sigColor,9}); text(65,36,'Slope','FontWeight','Bold');

changePosition(gca,[-0.03 0.09 0.09 -0.06]);
save_pdf(hf,mData.pdf_folder,'figure_scatter_centers_zMI.pdf',600);
return;
end

