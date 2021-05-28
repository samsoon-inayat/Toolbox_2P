function figure1_Distributions

% protocol_C = '10_C';
% protocol_A = '10_A';
ei_C = evalin('base','ei10_C');
ei_A = evalin('base','ei10_A');
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
% ET_C = evalin('base',sprintf('ET%s',protocol_C));
% ET_A = evalin('base',sprintf('ET%s',protocol_A));
selAnimals_C = 1:length(ei_C)
selAnimals_A = 1:length(ei_A)

% in the following variable all the measurements are in the matrices form
% for each variable colums indicate raster and stim marker types specified 
% the rows indicate condition numbers.
paramMs_C = parameter_matrices('get','10_C');
paramMs_A = parameter_matrices('get','10_A');
% after getting all matrics, we can apply selection criteria to select a
% subgroup of cells
% here is the selection criteria in make_selC_structure function
% cellsOrNot = 1; planeNumber = NaN; zMI_Th = 3; fwids = [1 120]; fcens = [0 140]; rs_th = 0.4;
cellsOrNot = 1; planeNumber = NaN; zMI_Th = NaN; fwids = NaN; fcens = NaN; rs_th = NaN; FR = NaN;
conditionsAndRasterTypes = [11 12];
selC = make_selC_struct(cellsOrNot,planeNumber,conditionsAndRasterTypes,zMI_Th,fwids,fcens,rs_th,NaN,NaN,FR);
[cpMs_C,pMs_C] = parameter_matrices('select','10_C',{paramMs_C,selC});
[cpMs_A,pMs_A] = parameter_matrices('select','10_A',{paramMs_A,selC});
% perc_cells_C = parameter_matrices('print','10_C',{cpMs_C,pMs_C,ET_C,selAnimals_C});
% perc_cells_A = parameter_matrices('print','10_A',{cpMs_A,pMs_A,ET_A,selAnimals_A});

all_conds = []; all_rts = [];
gAllVals_C = []; gAllVals_A = [];
for rr = 1:size(pMs_C,1)
    for cc = 1:size(pMs_C,2)
        tcond = conditionsAndRasterTypes(rr,cc);
        nds = dec2base(tcond,10) - '0';
        varNames{rr,cc} = sprintf('C%dR%d',nds(1),nds(2));
        all_conds = [all_conds nds(1)]; all_rts = [all_rts nds(2)];
        xticklabels{cc,rr} = sprintf('%s-%s',paramMs_C.stimMarkers{nds(2)},paramMs_C.rasterTypes{nds(2)}(1));
        for an = 1:length(selAnimals_C)
            zMIs_C(an,rr,cc) = nanmean(squeeze(pMs_C{rr,cc}.all_zMIs{selAnimals_C(an)}(nds(1),nds(2),:)));
            a_zMIs_C{an,rr,cc} = squeeze(pMs_C{rr,cc}.all_zMIs{selAnimals_C(an)}(nds(1),nds(2),:));
            gAllVals_C = [gAllVals_C;a_zMIs_C{an,rr,cc}];
        end
    end
end
for rr = 1:size(pMs_A,1)
    for cc = 1:size(pMs_A,2)
        tcond = conditionsAndRasterTypes(rr,cc);
        nds = dec2base(tcond,10) - '0';
        for an = 1:length(selAnimals_A)
            zMIs_A(an,rr,cc) = nanmean(squeeze(pMs_A{rr,cc}.all_zMIs{selAnimals_A(an)}(nds(1),nds(2),:)));
            a_zMIs_A{an,rr,cc} = squeeze(pMs_A{rr,cc}.all_zMIs{selAnimals_A(an)}(nds(1),nds(2),:));
            gAllVals_A = [gAllVals_A;a_zMIs_A{an,rr,cc}];
        end
    end
end
all_conds = unique(all_conds); all_rts = unique(all_rts);
% var_oi_A = squeeze(zMIs_A);
control = 1;
if control
zMIs = squeeze(a_zMIs_C);
% var_oi_C = squeeze(zMIs_C);
gAllVals = [gAllVals_C];
paramMs = paramMs_C;
colors = {'k',[0.5 0.5 0.5]};
else
    zMIs = squeeze(a_zMIs_A);
% var_oi_C = squeeze(zMIs_C);
gAllVals = [gAllVals_A];
paramMs = paramMs_A;
colors = {'r',[0.5 0 0]};
end
n = 0;
%%
runthis = 1;
if runthis
    data = zMIs;
    minBin = min(gAllVals);
    maxBin = 15;%max(gAllVals);
    incr = (maxBin-minBin)/100;
    
    hf = figure(1000);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 6 1.5 1],'color','w');
    hold on;
    [ha,hb,hca] = plotDistributions(data,'colors',colors,'maxY',90,'cumPos',[0.5 0.26 0.25 0.5],'min',minBin,'incr',incr,'max',maxBin);
    hold on;
    legs = [];
    for ii = 1:2%length(paramMs.stimMarkers)
        if ii == 1
            legs{ii} = sprintf('%s-%s (N = %d)',paramMs.stimMarkers{ii},paramMs.rasterTypes{ii}(1),size(data,1));
        else
            legs{ii} = sprintf('%s-%s',paramMs.stimMarkers{ii},paramMs.rasterTypes{ii}(1));
        end
    end
    ylim([0 100]);xlim([-3 8]);
    xlims = xlim; dx = xlims(2) - xlims(1); ylims = ylim; dy = ylims(2) - ylims(1);
    legs{ii+1} = [xlims(1)+dx/3 dx/30 ylims(1)+dy/3 dy/10];
    legs{1} = 'C1-AD (Air-Dist)'; legs{2} = 'C1-AT (Air-Time)'
    if size(data,2) > 2
        putLegend(ha,legs,'colors',colors,'sigR',{sigR,'ranova',sigColor,6});
    else
        putLegend(ha,legs,'colors',colors,'sigR',{[],'',sigColor,7});
    end
    axes(ha);
%     h = xlabel({'Mutual Information','Z-Score (zMI)'});%changePosition(h,[0 -dy/3 0]);
    h = xlabel('Mutual Information Z-Score (zMI)');changePosition(h,[-0.9 dy/25 0]);
    h = ylabel('Percentage');changePosition(h,[0.3 0 0]);
    set(gca,'FontSize',axes_font_size,'FontWeight','Bold');changePosition(ha,[0.07 0.13 -0.05 -0.15]);
    save_pdf(hf,mData.pdf_folder,sprintf('Distribution Of zMI %d _%d',all_conds(1),control),600);
   
end
