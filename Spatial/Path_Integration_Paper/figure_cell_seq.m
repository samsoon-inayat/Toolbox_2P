function figure_place_cells_vs_other_cells_1(fn,allRs,ccs)

ei = evalin('base','ei10');
mData = evalin('base','mData');

colors = mData.colors;
sigColor = mData.sigColor;
axes_font_size = mData.axes_font_size;

T = evalin('base','T10.T(selRecs,:)');

selAnimals = [1:9];
% in the following variable all the measurements are in the matrices form
% for each variable colums indicate raster and stim marker types specified 
% the rows indicate condition numbers.
paramMs = parameter_matrices('get');
% after getting all matrics, we can apply selection criteria to select a
% subgroup of cells
% here is the selection criteria in make_selC_structure function
% cellsOrNot = NaN; planeNumber = NaN; zMI_Th = 3; fwids = [0 140]; fcens = [0 140]; rs_th = 0.4;
cellsOrNot = NaN; planeNumber = NaN; zMI_Th = NaN; fwids = NaN; fcens = NaN; rs_th = NaN;
conditionsAndRasterTypes = [12 22 32 42]; selC = make_selC_struct(cellsOrNot,planeNumber,conditionsAndRasterTypes,zMI_Th,fwids,fcens,rs_th);
[cpMs,pMs] = parameter_matrices('select',{paramMs,selC});
% parameter_matrices('print percentages',{cpMs,pMs,T,selAnimals});

%%
all_trials = {1:2,3:4,5:7,8:10};%3:10;
all_trials = {1,2,3,4,5,6,7,8,9,10};%3:10;
sel_trials = [1 2 3 4];
% trials10 = 1;%3:9;
% align cells
stimMarkers = paramMs.stimMarkers;
rasterTypes = paramMs.rasterTypes;
CNi = 3; 
rasterTypeN = 1;

for ani = 1:length(selAnimals)
    an = selAnimals(ani);
    tei = ei(an);
    selCells = pMs{1}.cellSel{an};
    cns = paramMs.all_cns{an};
for crt = 1:length(conditionsAndRasterTypes)
    all_cellSeq = [];
for si = 1:length(all_trials)
    tcond = conditionsAndRasterTypes(crt);
    Ndigits = dec2base(tcond,10) - '0';
    mRsi = [];
%     for ani = 1:length(selAnimals)
% %         [si ani]
%         an = selAnimals(ani);
%         tei = ei(an);
%         selCells = pMs{1}.cellSel{an};
%         selCells = cpMs.cellSel{an};
%         cns = paramMs.all_cns{an};
        maxDistTime = paramMs.maxDistTime;
        if Ndigits(2) == 2
        [tempD cnso] = getParamValues('fromFrames.sp_rasters',tei,selC.plane_number,Ndigits(1),stimMarkers{Ndigits(2)},rasterTypes{Ndigits(2)},...
        cns(selCells,2:3),maxDistTime);
        else
        [tempD cnso] = getParamValues('rasters',tei,selC.plane_number,Ndigits(1),stimMarkers{Ndigits(2)},rasterTypes{Ndigits(2)},...
            cns(selCells,2:3),maxDistTime);
        end
        if length(tempD) == 0
            continue;
        end
        trials = all_trials{si};
         mR = findMeanRasters(tempD,trials);
        mRsi = [mRsi;mR];
%     end
    [temp,~,~] = getParamValues('',ei(1),1,1,stimMarkers{Ndigits(2)},rasterTypes{Ndigits(2)},'areCells',[Inf Inf]);
    dxs = diff(temp.xs); bin_width = dxs(1); xs = 0:bin_width:1000;
    allRs{si} = mRsi;
    time_xs{si} = xs(1:size(mRsi,2));
    raster_labels{si} = sprintf('Cond - %d, Rast - %d',Ndigits(1),Ndigits(2));
    [~,~,all_cellSeq(si,:)] = findPopulationVectorPlot(allRs{si},[]);
end

cell_pos = [];
for cc = 1:size(all_cellSeq,2)
    for rr = 1:size(all_cellSeq,1)
        cell_pos(cc,rr) = find(all_cellSeq(rr,:) == cc);
    end
end
% delta_cell_pos(crt,:) = mean(100*diff(cell_pos,1,2)/size(all_cellSeq,2),2);
% data{crt} = delta_cell_pos(crt,:);

cell_pos = 100*diff(cell_pos,1,2)/size(all_cellSeq,2);
mean_over_cells(crt,:) = mean(abs(cell_pos));

cell_pos = cell_pos - min(cell_pos(:));
cell_pos = cell_pos/max(cell_pos(:));
all_cell_pos{crt} = cell_pos;
e(crt) = entropy(cell_pos);
f(crt) = BoxCountfracDim(cell_pos);

end
all_mean_over_cells(:,:,ani) = mean_over_cells;
end
n = 0;
%%
figure(101);clf;hold on;
for ii = 1:crt
    plot(1:9,mean_over_cells(ii,:),'color',colors{ii});
%     subplot(1,crt,ii);
%     imagesc(all_cell_pos{ii});
end
return;
sigRM = significanceTesting(mean_over_cells');

n = 0;


%%
if 0
gAllVals = delta_cell_pos(:);
minBin = min(gAllVals);
maxBin = max(gAllVals);
incr = (maxBin-minBin)/20;
hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 6 4 3],'color','w');
hold on;
[ha,hb,hca,sigR] = plotDistributions(data,'colors',colors,'maxY',15,'cumPos',[0.5 0.26 0.25 0.5],...
    'min',minBin,'incr',incr,'max',maxBin,'do_mean','No');
hold on;
legs = [];
for ii = 1:length(conditionsAndRasterTypes)
    tcond = conditionsAndRasterTypes(ii);
    Ndigits = dec2base(tcond,10) - '0';
    legs{ii} = sprintf('%d-%d',Ndigits(1),Ndigits(2));
end
ylim([0 15]);
xlims = xlim; dx = xlims(2) - xlims(1); ylims = ylim; dy = ylims(2) - ylims(1);
legs{ii+1} = [xlims(1)+dx/4 dx/30 ylims(1)+dy/1 dy/15];
putLegend(ha,legs,'colors',colors,'sigR',{sigR,'anova',sigColor,10});
axes(ha);
h = xlabel('Mean Diff Pos');%changePosition(h,[0 -dy/3 0]);
h = ylabel('Percentage');changePosition(h,[-0 0 0]);
set(gca,'FontSize',axes_font_size+4,'FontWeight','Bold');changePosition(ha,[0.07 0.01 -0.05 -0.05]);
save_pdf(hf,mData.pdf_folder,sprintf('Distribution Of pos cell seq'),600);
end

%%
if 1
    
end
