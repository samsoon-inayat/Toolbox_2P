function motion_responsive_cells(ei)
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
if ~exist('ei','var')
    ei = evalin('base','d15_2(1)');
end

cai_sampling_rate = ei{1}.thorExp.frameRate;
effective_sampling_rate = 1/0.2;

pp = 1;

contexts = ei{1}.plane{pp}.contexts;
selContexts = [1 2 3 3 4 4 5 5 6 7];
samplingRate = {'Ca','Ef','Ef','Ef','Ef','Ef','Ef','Ef','Ca','Ef'};
rasterNames = {'light22T','air55T','air77T','airD','air77T','airD','air77T','airD','light22T','air55T'};
timeBefore = [2 4 4 NaN 4 NaN 4 NaN 2 4];

for ii = 1:length(selContexts)
    thisContext = contexts(selContexts(ii));
    disp(thisContext.name);
    if strcmp(samplingRate{ii},'Ca')
        cmdTxt = sprintf('tempR = thisContext.rasters.%s.fromFrames.sp_rasters;',rasterNames{ii});
        SR = cai_sampling_rate;
        cmdTxt1 = sprintf('tempDur = thisContext.rasters.%s.fromFrames.duration;',rasterNames{ii});
    else
        cmdTxt = sprintf('tempR = thisContext.rasters.%s.sp_rasters1;',rasterNames{ii});
        SR = effective_sampling_rate;
        cmdTxt1 = sprintf('tempDur = thisContext.rasters.%s.duration1;',rasterNames{ii});
    end
    eval(cmdTxt); eval(cmdTxt1);
    rasters{ii,1} = find_resp_mdata(tempR,tempDur,timeBefore(ii),SR,thisContext.name,ei);
end
n = 0;
%%
% plotRasters_multi(rasters,find(rasters{2}.resp.ps(:,1) < 0.05),[])
% Rs = rasters{6};plotRasters_simple(Rs,find(Rs.resp.p < 0.05));

meanRs = calc_mean_rasters(rasters,1:10);
meanRsTrials = calc_mean_rasters(rasters,{1:2,3:4,5:6,7:8,9:10});
meanRsRemap = calc_mean_rasters(rasters([2 10],1),1:10);
% meanRsRemap = calc_mean_rasters(rasters([1 9],1),1:10);
meanRsRemap = calc_mean_rasters(rasters([4 6 8],1),1:10);
%%
ind = 4;
ccs = rasters{ind}.resp.MIs > 3;
% ccs = rasters{2}.resp.p < 0.05;
sum(ccs)
% [popVs,pos_corr,cell_corr] = calc_pop_vector_corr(meanRs,ccs,[1,2]);
[popVsR,pos_corrR,cell_corrR] = calc_pop_vector_corr(meanRsRemap',ccs,[]);
[popVs,~,~] = calc_pop_vector_corr(meanRs,ccs,popVsR{1}.cell_nums');
sel_popVs = [1 2 4 6 8 9 10];
plot_pop_vectors(popVs(sel_popVs),rasters(sel_popVs));
n= 0 ;
%%
indr = 2;
Rsr = rasters{indr};
Rsp = rasters{2};
ccs = Rsr.resp.p < 0.05;
% ccs = Rsr.resp.MIs > 3;
sum(ccs)
[popVs,~,~] = calc_pop_vector_corr(meanRs,ccs,[indr,1]);
sel_popVs = [1 2 4 6 8 9 10];
plot_pop_vectors(popVs(sel_popVs),rasters(sel_popVs));

% plotRasters_simple(Rsp,find(Rsr.resp.p<0.05),[])

% plotRasters_simple(Rs,find(Rs.resp.MIs>3),[])
% plotRasters_multi(rasters,rastersD,find(resp.ps(:,1)<0.05),[])
% plotRasters_multi({rastersD1,rastersD2,rastersD3,rasters},find(zMIs1'>5),[])
%%
% sel_out = pos_corr{ind};
sel_out = pos_corrR{1};
ff = makeFigureRowsCols(107,[1 0.5 2 2],'RowsCols',size(sel_out),...
    'spaceRowsCols',[0.01 0.01],'rightUpShifts',[0.15 0.12],'widthHeightAdjustment',...
    [-85 -85]);
gg = 1;
set(gcf,'color','w');
set(gcf,'Position',[1 6 1.7 1.7]);
set(gcf,'Position',[1 2 6 6]);
FS = mData.axes_font_size;
maskdisp = triu(ones(5,5),0);

for rr = 1:size(sel_out,1)
    for cc = 1:size(sel_out,2)
        if ~maskdisp(rr,cc)
            delete(ff.h_axes(rr,cc));
            continue;
        end
        axes(ff.h_axes(rr,cc));
        this_pos_corr = sel_out{rr,cc};
        imagesc(this_pos_corr);
        box off;
        set(gca,'Ydir','Normal','linewidth',1,'FontSize',FS-1,'FontWeight','Bold');
        if rr == 1 && cc == 1
            cols = size(this_pos_corr,2);
            colsHalf = round(cols/2);
            ts = round(rasters{ind}.xs);
            set(gca,'YTick',[1 colsHalf cols],'YTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
            h = ylabel('Pos (cm)');%    changePosition(h,[0 0 0]);
            set(gca,'XTick',[]);
        elseif rr == size(sel_out,1) && cc == size(sel_out,2)
            h = xlabel('Pos (cm)');%    changePosition(h,[0 0 0]);
            cols = size(this_pos_corr,2);
            colsHalf = round(cols/2);
            ts = round(rasters{ind}.xs);
            set(gca,'XTick',[1 colsHalf cols],'XTickLabel',[ts(1) ts(colsHalf) ts(cols)]);
            set(gca,'YTick',[]);
        else
            axis off
        end
        minC = min(this_pos_corr(:));
        maxC = max(this_pos_corr(:));
%         text(5,cols-5,sprintf('(%.1f, %.1f)',minC,maxC),'FontSize',5,'Color','w');
%         hc = putColorBar(ff.h_axes(rr,cc),[0.0 0.03 0 -0.05],[minC maxC],5,'eastoutside',[0.07 0.11 0.1 0.16]);
    end
end
save_pdf(ff.hf,mData.pdf_folder,sprintf('pos corr trials'),600);
n = 0;
