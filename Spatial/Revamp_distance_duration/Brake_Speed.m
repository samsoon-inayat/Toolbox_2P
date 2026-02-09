magfac = mData.magfac; pdfopen = 0;
    ff = makeFigureRowsCols(108,[6 3 3.25 1.5],'RowsCols',[1 4],'spaceRowsCols',[0.03 -0.01],'rightUpShifts',[0.13 0.19],'widthHeightAdjustment',[10 -700]);
%     set(gcf,'color','w');    set(gcf,'Position',[10 3 3.5 3.75]);
    MY = 8; ysp = 1; mY = 0; % responsive cells
    stp = 0.22*magfac; widths = (0.5*ones(1,12)-0.08)*magfac; gap = 0.115*magfac; widths(3) = 0.84; widths(4) = 0.84;
    adjust_axes(ff,[mY MY],stp,widths,gap,{'Cell #'});
    % shift_axes(ff,[4 5 6;4 5 6;4 5 6],0.1,gap);
    % shift_axes_up(ff,[1 2 3 4;1 2 3 4],[0 0.15]);
    conf = repmat([2 7],1,2);

    o = o; G = 'C'; %all_cells_1 = all_cells(:,1:2:end); all_cells_1 = cell_list_op(all_cells_1,props.good_FR,'and');
    % all_cells_1 = cell_list_op(propsD.newMI.cells_R,props.good_FR,'and'); 
    pti = 4;
    paneltype = {'Time-Modulated Cells (PC)','Movement-Modulated Cells (PC)','Time-Modulated Cells (MI)','Movement-Modulated Cells (MI)'};
    poptype = {propsTB.newPC.cells_time,propsTB.newPC.cells_speed,propsTB.newMI.cells_time,propsTB.newMI.cells_speed};
    an = 4; all_cells_1 = poptype{pti}; %all_cells_1 = props.all;
    good_FR = all_cells_1(:,[1 3 2 4]);
    cbar_p_shift = [-0.011 0.09 -0.03 -0.23];
    tim = 1;
    si = [Ab_t_T Abs_t_T Ab_i_T Abs_i_T];
    Rs = o.Rs(:,si);mR = o.mR(:,si); props1 = get_props_Rs(Rs,ntrials);
    mRR = []; MRR = [];
    for ii = 1:size(ff.h_axes,2)
        tR = Rs{an,ii}; tmRR = tR.speed; thisRaster = tmRR;
        mRR(ii) = min(tmRR(:)); MRR(ii) = max(tmRR(:));
    end
    m_mRR = min(mRR); M_MRR = max(MRR);
    for ii = 1:size(ff.h_axes,2)
      % si = [Ab_t_T Abs_t_T Ab_i_T Abs_i_T];% 
    % Rs = o.Rs(:,si);mR = o.mR(:,si); props1 = get_props_Rs(Rs,ntrials);
        tR = Rs{an,ii}; 
        tmRR = tR.speed(:,1:(size(tR.speed,2)-1)); 
        % m_mRR = min(tmRR(:)); M_mRR = max(tmRR(:)); m_CRc = min(m_CRc_ii); M_CRc = max(M_CRc_ii); m_aCRc = min(m_aCRc_ii); M_aCRc = max(M_aCRc_ii);

        tR = Rs{an,ii};
        % NCells = size(tR.sp_rasters,3);
        thisRaster = tmRR;
        ax = ff.h_axes(1,ii);
        axes(ff.h_axes(1,ii));
        binwidth = tR.bin_width;
        totx = mxsz(ii)*binwidth;
        sz = size(thisRaster,2);
        xs2 = round(totx,0); xs1 = round(max(totx)/2,1);
        xdata = [0 xs1 xs2];
        ydata = [1 size(tmRR,1)];
        imagesc(xdata,ydata,tmRR,[m_mRR M_MRR]); set(gca,'Ydir','normal');
        if ii == 1
            ylabel('Trial #');
        end
        xlabel('Time (s)');
        
        set(gca,'YTick',[],'XTick',[0 xs2]);
        textstr = sprintf('C%d',conf(ii)); set_axes_top_text_no_line(ff.hf,gca,textstr,[0.0 0.065 0.05 0]);
        format_axes(gca);
        % mM = min(tmRR(:)); MM = max(tmRR(:)); 
        if ii == 4
          [hc,hca] = putColorBar(gca,cbar_p_shift,[m_mRR M_MRR],5,'eastoutside',[0.091 0.28 0.1 0.3]);
        end

        
         %******* make new axes and plot mean and gaussian fitting
        pos = get(ax,'Position'); ha = axes; set(ha,'units','inches');set(ha,'Position',pos);
        changePosition(ha,[0 pos(4)+0.15 0 -0.135]);
        has(ii) = ha;
    
        xs = linspace(0,sz*bw,sz);
        mSig = nanmean(tmRR);
        plot(xs,mSig,'b');hold on;
        box off;
        ylims = [-5 20];%ylim;
        set(gca,'XTick',[]);
        if ii == 1
            ylabel('Speed (cm/s)')
        end
        format_axes(gca);
    end
    colormap parula
    ht = set_axes_top_text(gcf,has(1,1:2),'AOn',{0.08,[0.1 -0.01 0 0]});set(ht,'FontWeight','Normal');
    ht = set_axes_top_text(gcf,has(1,3:4),'AOff',{0.08,[0.25 -0.01 0 0]});set(ht,'FontWeight','Normal');
    % ht = set_axes_top_text_no_line(gcf,ff.h_axes(1,3:4),paneltype{pti},[-0.15 0.18 0 0]);set(ht,'FontWeight','Normal','FontSize',7.25);
    filename = save_pdf(ff.hf,mData.pdf_folder,sprintf('Brake_speed_%d_%d.pdf',ntrials,tim),600);
if pdfopen
winopen(filename)
end

%%
si = [Ab_t_T Abs_t_T Ar_t_T ArL_t_T Ars_t_T Ab_i_T Abs_i_T Ar_i_T ArL_i_T Ars_i_T];
Rs = o.Rs(:,si);mR = o.mR(:,si); props1 = get_props_Rs(Rs,ntrials);

mspeed = [];
for an = 1:5
    for ii = 1:length(si)
        tR = Rs{an,ii};
        sRaster = tR.speed;
        mspeed(an,ii) = nanmean(nanmean(sRaster,2));
    end
end

avar = mspeed;

fac_names = {'AP','CN'}; fac_levels = [2,5];
% fac_names = {'CN','AP','BT','RF','CT'}; fac_levels = [3,2,2,3,4];
% fac_names = {'CN','AP','BT','CT'}; fac_levels = [3,2,2,7];
[within,dvn,xlabels,awithinD] = make_within_table(fac_names,fac_levels);
dataT = make_between_table({avar},dvn);
ra = RMA(dataT,within,{0.05,{''}});
clc
print_for_manuscript(ra)

%% visualizing the results in the previous section
tcolors = repmat(mData.dcolors(1:5),1,2); MY = 50; ysp = 1; mY = 0; ystf = 2; ysigf = 0.025;titletxt = ''; ylabeltxt = {'Cells (%)'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova([],ra,{'AP:CN','hsd',0.05},[1 1.75],tcolors,[mY MY ysp ystf ysigf],mData);
%% Figure 5G ... comparison of Brake speed vs no-brake speed
% visualization
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size; dcolors = mData.dcolors;
tcolors = repmat(mData.dcolors(1:5),1,2);
% figure(300);clf; ha = gca;
ff = makeFigureRowsCols(2020,[5 4 3.5 1.5],'RowsCols',[1 1],'spaceRowsCols',[0.07 0],'rightUpShifts',[0.1 0.25],'widthHeightAdjustment',[-120 -300]);
MY = 50; ysp = 2.75; mY = 0; ystf = 5; ysigf = 0.05;titletxt = ''; ylabeltxt = {'PDF'}; % for all cells (vals) MY = 80
[hbs,xdata,mVar,semVar,combs,p,h] = view_results_rmanova(ff.h_axes(1,1),ra,{'AP:CN','hsd',0.05},[1 1.5],tcolors,[mY MY ysp ystf ysigf],mData);
% make_bars_hollow(hbs(2))
set(gca,'xcolor','k','ycolor','k','xlim',xlim,'ylim',ylim,...
    'XTick',xdata,'XTickLabel',{'C2','C7','C3','C4','C5'});xtickangle(20);
ylabel({'Avg. Speed (cm/s)'});
set_bar_graph_sub_xtick_text(ff.hf,gca,hbs,5,{'AOn','Air-Off'},{[0 0]});
format_axes(gca);
% ht = set_axes_top_text_no_line(ff.hf,gca,sprintf('Pooled'),[0.051 0.0 0 0]); 
save_pdf(ff.hf,mData.pdf_folder,sprintf('bar_graphs.pdf'),600);