%% population vector and correlation sensory
    magfac = mData.magfac;
    ff = makeFigureRowsCols(108,[6 3 2.6 2.5],'RowsCols',[3 3],'spaceRowsCols',[0.03 -0.01],'rightUpShifts',[0.1 0.11],'widthHeightAdjustment',[10 -85]);
%     set(gcf,'color','w');    set(gcf,'Position',[10 3 3.5 3.75]);
    MY = 8; ysp = 1; mY = 0; % responsive cells
    stp = 0.15*magfac; widths = ([0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5]+0.11)*magfac; gap = 0.19*magfac;
    adjust_axes(ff,[mY MY],stp,widths,gap,{'Cell #'});
    
    an = 4;o = o; G = 'C';  
%     an = 4;o = oC; G = 'C';
%     si = [C1_t_D C2_t_D C3_t_D C4_t_D];
    si = [Ar_t_D ArL_t_D Ars_t_D];% 
%     si = [Ar_B_T ArL_B_T Ars_B_T];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    good_FR = cell_list_op(props1,{'vals','good_FR'}); 
    good_FR = cell_list_op(props1,{'vals'}); 
% %     good_FR = dzMI_FD.resp_T_g_D;
%     good_FR = cell_list_op(props1,{'vals'});
    good_FR = dis_cells_T;
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,good_FR,0);
    mRRm = [];
    for ii = 1:size(mRR,2)
        m_mRR_ii(ii) = min(min(mRR{an,ii}));        M_mRR_ii(ii) = max(max(mRR{an,ii}));
        m_CRc_ii(ii) = min(min(CRc{an,ii}));        M_CRc_ii(ii) = max(max(CRc{an,ii}));
        m_aCRc_ii(ii) = min(min(aCRc{ii}));        M_aCRc_ii(ii) = max(max(aCRc{ii}));
    end
    m_mRR = min(m_mRR_ii); M_mRR = max(M_mRR_ii); m_CRc = min(m_CRc_ii); M_CRc = max(M_CRc_ii); m_aCRc = min(m_aCRc_ii); M_aCRc = max(M_aCRc_ii);
    cbar_p_shift = [-0.011 0.09 -0.03 -0.3];
    for ii = 1:size(ff.h_axes,2)
        tR = Rs{an,ii};
        NCells = size(tR.sp_rasters,3);
        tmRR = mRR{an,ii};
        axes(ff.h_axes(1,ii));
        if ii < 5
            xdata = [0 75 150];
        else
            xdata = [0 7.5 15];
        end
        ydata = [1 size(tmRR,1)];
        imagesc(xdata,ydata,tmRR,[m_mRR M_mRR]); set(gca,'Ydir','normal');
        if ii == 1
            ylabel('Cell #');
        end
        set(gca,'YTick',[],'XTick',[]);

        ylims = ylim;
        if ii == 1
            textstr = sprintf('%d/%d',floor(ylims(2)),NCells); set_axes_top_text_no_line(ff.hf,gca,textstr,[0 -0.125 0 0]);
        else
            textstr = sprintf('%d',floor(ylims(2))); set_axes_top_text_no_line(ff.hf,gca,textstr,[0 -0.125 0 0]);
        end
        textstr = sprintf('Config. %d (C%d)',ii+2,ii+2); set_axes_top_text_no_line(ff.hf,gca,textstr,[0 -0.093 0 0]);
        format_axes(gca);
        mM = min(tmRR(:)); MM = max(tmRR(:)); [hc,hca] = putColorBar(gca,cbar_p_shift,[mM MM],5,'eastoutside',[0.09 0.11 0.09 0.16]);
        
        axes(ff.h_axes(2,ii));
        tmRR = CRc{an,ii};
        if sum(isnan(tmRR(:))) > 0
            tmRR = fillmissing(tmRR,'linear',2,'EndValues','nearest');
            if sum(isnan(tmRR(:))) > 0
                tmRR = tmRR'; tmRR = fillmissing(tmRR,'linear',2,'EndValues','nearest'); tmRR = tmRR';
            end
        end
        
        if ii < 5
            xdata = [0 75 150];
        else
            xdata = [0 7.5 15];
        end
        ydata = [1 size(tmRR,1)];
        imagesc(xdata,ydata,tmRR,[m_CRc M_CRc]); set(gca,'Ydir','normal');
        if ii == 1
            ylabel('Distance (cm)');
        end
        set(gca,'YTick',[],'XTick',[]);

        ylims = ylim;
        textstr = sprintf('%d',ceil(ylims(2)));
%         set_axes_top_text_no_line(ff.hf,gca,textstr,[0 -0.07 0 0]);
        format_axes(gca);
        mM = -0.3;min(tmRR(:)); MM = max(tmRR(:)); [hc,hca] = putColorBar(gca,cbar_p_shift,[mM MM],5,'eastoutside',[0.09 0.11 0.09 0.16]);
        
        axes(ff.h_axes(3,ii));
        tmRR = aCRc{ii};
        if ii < 5
            xdata = [0 75 150];
        else
            xdata = [0 7.5 15];
        end
        ydata = [1 size(tmRR,1)];
        imagesc(xdata,ydata,tmRR,[m_aCRc M_aCRc]); set(gca,'Ydir','normal');
        if ii == 1
            ylabel('Distance (cm)');
        end
        set(gca,'YTick',[]);
        if ii < 5
            xlabel('Distance (cm)');
        else
            xlabel('Time (sec)');
        end
        ylims = ylim;
        textstr = sprintf('%d',ceil(ylims(2)));
%         set_axes_top_text_no_line(ff.hf,gca,textstr,[0 -0.07 0 0]);
        format_axes(gca);
        mM = -0.3;min(tmRR(:)); MM = max(tmRR(:)); [hc,hca] = putColorBar(gca,cbar_p_shift,[mM MM],5,'eastoutside',[0.09 0.11 0.09 0.16]);
    end
    colormap jet
    save_pdf(ff.hf,mData.pdf_folder,sprintf('pop_vectors_%d_%c.pdf',ntrials,G),600);

%% population vector and correlation sensory for time rasters

    magfac = mData.magfac;
    ff = makeFigureRowsCols(108,[6 3 2.6 2.5],'RowsCols',[3 3],'spaceRowsCols',[0.03 -0.01],'rightUpShifts',[0.1 0.11],'widthHeightAdjustment',[10 -85]);
%     set(gcf,'color','w');    set(gcf,'Position',[10 3 3.5 3.75]);
    MY = 8; ysp = 1; mY = 0; % responsive cells
    stp = 0.15*magfac; widths = ([0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5]+0.11)*magfac; gap = 0.19*magfac;
    adjust_axes(ff,[mY MY],stp,widths,gap,{'Cell #'});
    
    an = 4;o = o; G = 'C';  
    si = [Ar_i_T ArL_i_T Ars_i_T];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 60;
    props1 = get_props_Rs(Rs,ntrials);
    good_FR = cell_list_op(props1,{'good_zMI','good_Gauss'});
    good_FR = cell_list_op(props1,{'vals','good_FR'});
    good_FR = dur_cells_I;
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,good_FR,0);
    mRRm = [];
    for ii = 1:size(mRR,2)
        m_mRR_ii(ii) = min(min(mRR{an,ii}));        M_mRR_ii(ii) = max(max(mRR{an,ii}));
        m_CRc_ii(ii) = min(min(CRc{an,ii}));        M_CRc_ii(ii) = max(max(CRc{an,ii}));
        m_aCRc_ii(ii) = min(min(aCRc{ii}));        M_aCRc_ii(ii) = max(max(aCRc{ii}));
    end
    m_mRR = min(m_mRR_ii); M_mRR = max(M_mRR_ii); m_CRc = min(m_CRc_ii); M_CRc = max(M_CRc_ii); m_aCRc = min(m_aCRc_ii); M_aCRc = max(M_aCRc_ii);
    cbar_p_shift = [-0.011 0.09 -0.03 -0.3];
    for ii = 1:size(ff.h_axes,2)
        tR = Rs{an,ii};
        NCells = size(tR.sp_rasters,3);
        tmRR = mRR{an,ii};
        axes(ff.h_axes(1,ii));
        if ii < 5
            xdata = [0 75 150];
        else
            xdata = [0 7.5 15];
        end
        ydata = [1 size(tmRR,1)];
        imagesc(xdata,ydata,tmRR,[m_mRR M_mRR]); set(gca,'Ydir','normal');
        if ii == 1
            ylabel('Cell #');
        end
        set(gca,'YTick',[],'XTick',[]);

        ylims = ylim;
        if ii == 1
            textstr = sprintf('%d/%d',floor(ylims(2)),NCells); set_axes_top_text_no_line(ff.hf,gca,textstr,[0 -0.125 0 0]);
        else
            textstr = sprintf('%d',floor(ylims(2))); set_axes_top_text_no_line(ff.hf,gca,textstr,[0 -0.125 0 0]);
        end
%         textstr = sprintf('%d',floor(ylims(2)));set_axes_top_text_no_line(ff.hf,gca,textstr,[0 -0.13 0 0]);
        textstr = sprintf('Config. %d (C%d)',ii+2,ii+2); set_axes_top_text_no_line(ff.hf,gca,textstr,[0 -0.093 0 0]);
        format_axes(gca);
        mM = min(tmRR(:)); MM = max(tmRR(:)); [hc,hca] = putColorBar(gca,cbar_p_shift,[mM MM],5,'eastoutside',[0.09 0.11 0.09 0.16]);
        
        axes(ff.h_axes(2,ii));
        tmRR = CRc{an,ii};
        if sum(isnan(tmRR(:))) > 0
            tmRR = fillmissing(tmRR,'linear',2,'EndValues','nearest');
            if sum(isnan(tmRR(:))) > 0
                tmRR = tmRR'; tmRR = fillmissing(tmRR,'linear',2,'EndValues','nearest'); tmRR = tmRR';
            end
        end
        
            xdata = [0 7.5 15];

        ydata = [1 size(tmRR,1)];
        imagesc(xdata,ydata,tmRR,[m_CRc M_CRc]); set(gca,'Ydir','normal');
        if ii == 1
            ylabel('Time (s)');
        end
        set(gca,'YTick',[],'XTick',[]);

        ylims = ylim;
        textstr = sprintf('%d',ceil(ylims(2)));
%         set_axes_top_text_no_line(ff.hf,gca,textstr,[0 -0.07 0 0]);
        format_axes(gca);
        mM = -0.3;min(tmRR(:)); MM = max(tmRR(:)); [hc,hca] = putColorBar(gca,cbar_p_shift,[mM MM],5,'eastoutside',[0.09 0.11 0.09 0.16]);
        
        axes(ff.h_axes(3,ii));
        tmRR = aCRc{ii};

            xdata = [0 7.5 15];

        ydata = [1 size(tmRR,1)];
        imagesc(xdata,ydata,tmRR,[m_aCRc M_aCRc]); set(gca,'Ydir','normal');
        if ii == 1
            ylabel('Time (s)');
        end
        set(gca,'YTick',[]);
 
            xlabel('Time (s)');
   
        ylims = ylim;
        textstr = sprintf('%d',ceil(ylims(2)));
%         set_axes_top_text_no_line(ff.hf,gca,textstr,[0 -0.07 0 0]);
        format_axes(gca);
        mM = -0.3;min(tmRR(:)); MM = max(tmRR(:)); [hc,hca] = putColorBar(gca,cbar_p_shift,[mM MM],5,'eastoutside',[0.09 0.11 0.09 0.16]);
    end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('pop_vectors_%d_%c.pdf',ntrials,G),600);
   
%% population vector and correlation sensory (Dist - time) change var tim
    magfac = mData.magfac; pdfopen = 1;
    ff = makeFigureRowsCols(108,[6 3 3.5 2.15],'RowsCols',[3 6],'spaceRowsCols',[0.03 -0.01],'rightUpShifts',[0.1 0.12],'widthHeightAdjustment',[10 -140]);
%     set(gcf,'color','w');    set(gcf,'Position',[10 3 3.5 3.75]);
    MY = 8; ysp = 1; mY = 0; % responsive cells
    stp = 0.15*magfac; widths = (0.5*ones(1,12)-0.08)*magfac; gap = 0.115*magfac;
    adjust_axes(ff,[mY MY],stp,widths,gap,{'Cell #'});
    shift_axes(ff,[4 5 6;4 5 6;4 5 6],0.1,gap);
    shift_axes_up(ff,[1 2 3 4 5 6;1 2 3 4 5 6],[0 0.15 0 0]);
    conf = repmat([3 4 5],1,2);
    o = o; G = 'C'; %all_cells_1 = all_cells(:,1:2:end); all_cells_1 = cell_list_op(all_cells_1,props.good_FR,'and');
    % all_cells_1 = cell_list_op(propsD.newMI.cells_R,props.good_FR,'and'); 
    an = 4; all_cells_1 = propsT.newMI.cells_speed; %all_cells_1 = props.all;
    good_FR = all_cells_1(:,[1 3 5 2 4 6]);
    cbar_p_shift = [-0.011 0.09 -0.03 -0.3];
    for ii = 1:size(ff.h_axes,2)
            if ii < 4
                tim = 0;
            else
                tim = 1;
            end
                    if tim
                      si = [Ar_t_T ArL_t_T Ars_t_T Ar_i_T ArL_i_T Ars_i_T];% 
                    else
                      si = [Ar_t_D ArL_t_D Ars_t_D Ar_i_D ArL_i_D Ars_i_D];% 
                    end
                    Rs = o.Rs(:,si);mR = o.mR(:,si); props1 = get_props_Rs(Rs,ntrials);
                    [CRc,aCRc,mRR,~,mxsz] = find_population_vector_corr(Rs,mR,good_FR,0);
                    mRRm = [];
                    for iii = 1:size(mRR,2)
                        m_mRR_ii(iii) = min(min(mRR{an,iii}));        M_mRR_ii(iii) = max(max(mRR{an,iii}));
                        m_CRc_ii(iii) = min(min(CRc{an,iii}));        M_CRc_ii(iii) = max(max(CRc{an,iii}));
                        m_aCRc_ii(iii) = min(min(aCRc{iii}));        M_aCRc_ii(iii) = max(max(aCRc{iii}));
                    end
                    m_mRR = min(m_mRR_ii); M_mRR = max(M_mRR_ii); m_CRc = min(m_CRc_ii); M_CRc = max(M_CRc_ii); m_aCRc = min(m_aCRc_ii); M_aCRc = max(M_aCRc_ii);

        tR = Rs{an,ii};
        NCells = size(tR.sp_rasters,3);
        tmRR = mRR{an,ii};
        axes(ff.h_axes(1,ii));
        binwidth = tR.bin_width;
        totx = mxsz(ii)*binwidth;
        xs2 = round(totx,0); xs1 = round(max(totx)/2,1);
        xdata = [0 xs1 xs2];
        ydata = [1 size(tmRR,1)];
        imagesc(xdata,ydata,tmRR,[m_mRR M_mRR]); set(gca,'Ydir','normal');
        if ii == 1
            ylabel('Cell #');
        end
        set(gca,'YTick',[],'XTick',[]);

        ylims = ylim;
        if ii == 1
            textstr = sprintf('C%d  %d/%d',conf(ii),floor(ylims(2)),NCells); set_axes_top_text_no_line(ff.hf,gca,textstr,[0.0 0.051 0.05 0]);
        else
            if sum(tmRR,'all') == 0
                textstr = sprintf('C%d    %d',conf(ii),0); set_axes_top_text_no_line(ff.hf,gca,textstr,[0.0 0.051 0 0]);
            else
                textstr = sprintf('C%d    %d',conf(ii),floor(ylims(2))); set_axes_top_text_no_line(ff.hf,gca,textstr,[0.0 0.051 0 0]);
            end
        end
%         textstr = sprintf('C%d',conf(ii)); set_axes_top_text_no_line(ff.hf,gca,textstr,[0 -0.085 0 0]);
        format_axes(gca);
        if ii == 6
          mM = min(tmRR(:)); MM = max(tmRR(:)); 
          [hc,hca] = putColorBar(gca,cbar_p_shift,[m_mRR M_mRR],5,'eastoutside',[0.15 0.18 0.12 0.3]);
        end
        
        axes(ff.h_axes(2,ii));
        tmRR = CRc{an,ii};
        if sum(isnan(tmRR(:))) > 0
            tmRR = fillmissing(tmRR,'linear',2,'EndValues','nearest');
            if sum(isnan(tmRR(:))) > 0
                tmRR = tmRR'; tmRR = fillmissing(tmRR,'linear',2,'EndValues','nearest'); tmRR = tmRR';
            end
        end
        if tim == 0
          if ii <= 3
              totx = 150;
          else
          totx = max(tR.xs);%mxsz(ii)*binwidth;
          end
        else
          if ii > 3
              totx = 15;
          else
          totx = max(tR.xs);%mxsz(ii)*binwidth;
          end
        end
        xs2 = round(totx,0); xs1 = round(max(totx)/2,1);
        xdata = [0 xs2];
        ydata = [1 size(tmRR,1)];
        imagesc(xdata,ydata,tmRR,[m_CRc M_CRc]); set(gca,'Ydir','normal');
        if ii == 1
          if tim
            ylabel('Time (s)');
          else
            ylabel('Dist (cm)');
          end
        end
        set(gca,'YTick',[],'XTick',xdata);

        ylims = ylim;
        textstr = sprintf('%d',ceil(ylims(2)));
%         set_axes_top_text_no_line(ff.hf,gca,textstr,[0 -0.07 0 0]);
        format_axes(gca);
        mM = -0.3;min(tmRR(:)); MM = max(tmRR(:));
        if ii == 6
          [hc,hca] = putColorBar(gca,cbar_p_shift,[m_CRc M_CRc],5,'eastoutside',[0.13 0.18 0.12 0.3]);
        end
        
        axes(ff.h_axes(3,ii));
        tmRR = aCRc{ii};
        if tim
          if ii > 3
              xs1 = 7.5; xs2 = 15;
          else
              totx = mxsz(ii)*binwidth;
              xs2 = round(totx,0); xs1 = round(max(totx)/2,1);
          end
        else
          if ii < 4
              xs1 = 75; xs2 = 150;
          else
              totx = mxsz(ii)*binwidth;
              xs2 = round(totx,0); xs1 = round(max(totx)/2,1);
          end
        end
        xdata = [0 xs1 xs2];
        ydata = [1 size(tmRR,1)];
        imagesc(xdata,ydata,tmRR,[m_aCRc M_aCRc]); set(gca,'Ydir','normal');
        if ii == 1
            if tim
            ylabel('Time (s)');
          else
            ylabel('Dist (cm)');
          end
        end
        set(gca,'YTick',[],'XTick',[0 xs2]);
        if tim
            xlabel('Time (s)');
          else
            xlabel('Dist (cm)');
          end
        
        ylims = ylim;
        textstr = sprintf('%d',ceil(ylims(2)));
%         set_axes_top_text_no_line(ff.hf,gca,textstr,[0 -0.07 0 0]);
        format_axes(gca);
        mM = -0.3;min(tmRR(:)); MM = max(tmRR(:)); 
        if ii == 6
          [hc,hca] = putColorBar(gca,cbar_p_shift,[m_aCRc M_aCRc],5,'eastoutside',[0.13 0.18 0.12 0.3]);
        end
    end
    colormap jet
    ht = set_axes_top_text(gcf,ff.h_axes(1,1:3),'AOn',{0.08,[0.2 -0.05 0 0]});set(ht,'FontWeight','Bold');
    ht = set_axes_top_text(gcf,ff.h_axes(1,4:6),'AOff',{0.08,[0.17 -0.05 0 0]});set(ht,'FontWeight','Bold');
    filename = save_pdf(ff.hf,mData.pdf_folder,sprintf('pop_vectors_%d_%d.pdf',ntrials,tim),600);
if pdfopen
winopen(filename)
end

%% population vector and correlation sensory (Dist - time) change var tim
    magfac = mData.magfac;
    ff = makeFigureRowsCols(108,[6 3 3.5 2.15],'RowsCols',[3 6],'spaceRowsCols',[0.03 -0.01],'rightUpShifts',[0.1 0.12],'widthHeightAdjustment',[10 -140]);
%     set(gcf,'color','w');    set(gcf,'Position',[10 3 3.5 3.75]);
    MY = 8; ysp = 1; mY = 0; % responsive cells
    stp = 0.15*magfac; widths = (0.5*ones(1,12)-0.08)*magfac; gap = 0.115*magfac;
    adjust_axes(ff,[mY MY],stp,widths,gap,{'Cell #'});
    shift_axes(ff,[4 5 6;4 5 6;4 5 6],0.1,gap);
    shift_axes_up(ff,[1 2 3 4 5 6;1 2 3 4 5 6],[0 0.15 0 0]);
    conf = repmat([3 4 5],1,2);
    o = o; G = 'C';  
    an = 3; tim = 1; good_FR = good_FR_t;
    if tim
      si = [Ar_t_T ArL_t_T Ars_t_T Ar_i_T ArL_i_T Ars_i_T];% 
      Rs = o.Rs(:,si);mR = o.mR(:,si); props1 = get_props_Rs(Rs,ntrials);
      % good_FR = [dur_cells_T dur_cells_I];
%           good_FR = cell_list_op(props1,{'vals'}); 
    else
      si = [Ar_t_D ArL_t_D Ars_t_D Ar_i_D ArL_i_D Ars_i_D];% 
      Rs = o.Rs(:,si);mR = o.mR(:,si); props1 = get_props_Rs(Rs,ntrials);
      % good_FR = [dis_cells_T dis_cells_I];
%           good_FR = cell_list_op(props1,{'vals'}); 
    end
    [CRc,aCRc,mRR,~,mxsz] = find_population_vector_corr(Rs,mR,good_FR,0);
    mRRm = [];
    for ii = 1:size(mRR,2)
        m_mRR_ii(ii) = min(min(mRR{an,ii}));        M_mRR_ii(ii) = max(max(mRR{an,ii}));
        m_CRc_ii(ii) = min(min(CRc{an,ii}));        M_CRc_ii(ii) = max(max(CRc{an,ii}));
        m_aCRc_ii(ii) = min(min(aCRc{ii}));        M_aCRc_ii(ii) = max(max(aCRc{ii}));
    end
    m_mRR = min(m_mRR_ii); M_mRR = max(M_mRR_ii); m_CRc = min(m_CRc_ii); M_CRc = max(M_CRc_ii); m_aCRc = min(m_aCRc_ii); M_aCRc = max(M_aCRc_ii);
    cbar_p_shift = [-0.011 0.09 -0.03 -0.3];
    for ii = 1:size(ff.h_axes,2)
        tR = Rs{an,ii};
        NCells = size(tR.sp_rasters,3);
        tmRR = mRR{an,ii};
        axes(ff.h_axes(1,ii));
        binwidth = tR.bin_width;
        totx = mxsz(ii)*binwidth;
        xs2 = round(totx,0); xs1 = round(max(totx)/2,1);
        xdata = [0 xs1 xs2];
        ydata = [1 size(tmRR,1)];
        imagesc(xdata,ydata,tmRR,[m_mRR M_mRR]); set(gca,'Ydir','normal');
        if ii == 1
            ylabel('Cell #');
        end
        set(gca,'YTick',[],'XTick',[]);

        ylims = ylim;
        if ii == 1
            textstr = sprintf('C%d  %d/%d',conf(ii),floor(ylims(2)),NCells); set_axes_top_text_no_line(ff.hf,gca,textstr,[0.0 0.051 0.05 0]);
        else
            textstr = sprintf('C%d    %d',conf(ii),floor(ylims(2))); set_axes_top_text_no_line(ff.hf,gca,textstr,[0.0 0.051 0 0]);
        end
%         textstr = sprintf('C%d',conf(ii)); set_axes_top_text_no_line(ff.hf,gca,textstr,[0 -0.085 0 0]);
        format_axes(gca);
        if ii == 6
          mM = min(tmRR(:)); MM = max(tmRR(:)); 
          [hc,hca] = putColorBar(gca,cbar_p_shift,[m_mRR M_mRR],5,'eastoutside',[0.15 0.18 0.12 0.3]);
        end
        
        axes(ff.h_axes(2,ii));
        tmRR = CRc{an,ii};
        if sum(isnan(tmRR(:))) > 0
            tmRR = fillmissing(tmRR,'linear',2,'EndValues','nearest');
            if sum(isnan(tmRR(:))) > 0
                tmRR = tmRR'; tmRR = fillmissing(tmRR,'linear',2,'EndValues','nearest'); tmRR = tmRR';
            end
        end
        if tim == 0
          if ii <= 3
              totx = 150;
          else
          totx = max(tR.xs);%mxsz(ii)*binwidth;
          end
        else
          if ii > 3
              totx = 15;
          else
          totx = max(tR.xs);%mxsz(ii)*binwidth;
          end
        end
        xs2 = round(totx,0); xs1 = round(max(totx)/2,1);
        xdata = [0 xs2];
        ydata = [1 size(tmRR,1)];
        imagesc(xdata,ydata,tmRR,[m_CRc M_CRc]); set(gca,'Ydir','normal');
        if ii == 1
          if tim
            ylabel('Time (s)');
          else
            ylabel('Dist (cm)');
          end
        end
        set(gca,'YTick',[],'XTick',xdata);

        ylims = ylim;
        textstr = sprintf('%d',ceil(ylims(2)));
%         set_axes_top_text_no_line(ff.hf,gca,textstr,[0 -0.07 0 0]);
        format_axes(gca);
        mM = -0.3;min(tmRR(:)); MM = max(tmRR(:));
        if ii == 6
          [hc,hca] = putColorBar(gca,cbar_p_shift,[m_CRc M_CRc],5,'eastoutside',[0.13 0.18 0.12 0.3]);
        end
        
        axes(ff.h_axes(3,ii));
        tmRR = aCRc{ii};
        if tim
          if ii > 3
              xs1 = 7.5; xs2 = 15;
          else
              totx = mxsz(ii)*binwidth;
              xs2 = round(totx,0); xs1 = round(max(totx)/2,1);
          end
        else
          if ii < 4
              xs1 = 75; xs2 = 150;
          else
              totx = mxsz(ii)*binwidth;
              xs2 = round(totx,0); xs1 = round(max(totx)/2,1);
          end
        end
        xdata = [0 xs1 xs2];
        ydata = [1 size(tmRR,1)];
        imagesc(xdata,ydata,tmRR,[m_aCRc M_aCRc]); set(gca,'Ydir','normal');
        if ii == 1
            if tim
            ylabel('Time (s)');
          else
            ylabel('Dist (cm)');
          end
        end
        set(gca,'YTick',[],'XTick',[0 xs2]);
        if tim
            xlabel('Time (s)');
          else
            xlabel('Dist (cm)');
          end
        
        ylims = ylim;
        textstr = sprintf('%d',ceil(ylims(2)));
%         set_axes_top_text_no_line(ff.hf,gca,textstr,[0 -0.07 0 0]);
        format_axes(gca);
        mM = -0.3;min(tmRR(:)); MM = max(tmRR(:)); 
        if ii == 6
          [hc,hca] = putColorBar(gca,cbar_p_shift,[m_aCRc M_aCRc],5,'eastoutside',[0.13 0.18 0.12 0.3]);
        end
    end
    colormap jet
    ht = set_axes_top_text(gcf,ff.h_axes(1,1:3),'AOn',{0.08,[0.2 -0.05 0 0]});set(ht,'FontWeight','Bold');
    ht = set_axes_top_text(gcf,ff.h_axes(1,4:6),'AOff',{0.08,[0.17 -0.05 0 0]});set(ht,'FontWeight','Bold');
    save_pdf(ff.hf,mData.pdf_folder,sprintf('pop_vectors_%d_%d.pdf',ntrials,tim),600);

%% population vector and correlation sensory (Dist and time) conjunctive
    magfac = mData.magfac;
    ff = makeFigureRowsCols(108,[6 3 3.5 2.15],'RowsCols',[3 4],'spaceRowsCols',[0.03 -0.01],'rightUpShifts',[0.1 0.12],'widthHeightAdjustment',[10 -140]);
%     set(gcf,'color','w');    set(gcf,'Position',[10 3 3.5 3.75]);
    MY = 8; ysp = 1; mY = 0; % responsive cells
    stp = 0.15*magfac; widths = (0.5*ones(1,12)-0.08)*magfac; gap = 0.115*magfac;
    adjust_axes(ff,[mY MY],stp,widths,gap,{'Cell #'});
%     shift_axes(ff,[4 5 6;4 5 6;4 5 6],0.1,gap);
%     shift_axes_up(ff,[1 2 3 4 5 6;1 2 3 4 5 6],[0 0.15 0 0]);
    conf = repmat([3 4 5],1,2);
    an = 3; o = o; G = 'C';  
    cni = 1;

    si = [Ar_t_T ArL_t_T Ars_t_T Ar_i_T ArL_i_T Ars_i_T Ar_t_D ArL_t_D Ars_t_D Ar_i_D ArL_i_D Ars_i_D];% 
    si = si([([1 4]+(cni-1)) ([1 4]+(cni-1))+6]);
%     good_FR = [dur_cells_T(:,cni) dur_cells_I(:,cni) dis_cells_T(:,cni) dis_cells_I(:,cni)];

    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    good_FR = cell_list_op(props1,{'vals'}); 
%     good_FR = cell_list_op(dis_cells_T,dis_cells_I,'and',1);
    good_FR = repmat(good_FR,1,6);
    
    [CRc,aCRc,mRR,~,mxsz] = find_population_vector_corr(Rs,mR,good_FR,1);
    mRRm = [];
    for ii = 1:size(mRR,2)
        m_mRR_ii(ii) = min(min(mRR{an,ii}));        M_mRR_ii(ii) = max(max(mRR{an,ii}));
        m_CRc_ii(ii) = min(min(CRc{an,ii}));        M_CRc_ii(ii) = max(max(CRc{an,ii}));
        m_aCRc_ii(ii) = min(min(aCRc{ii}));        M_aCRc_ii(ii) = max(max(aCRc{ii}));
    end
    m_mRR = min(m_mRR_ii); M_mRR = max(M_mRR_ii); m_CRc = min(m_CRc_ii); M_CRc = max(M_CRc_ii); m_aCRc = min(m_aCRc_ii); M_aCRc = max(M_aCRc_ii);
    cbar_p_shift = [-0.011 0.09 -0.03 -0.3];
    for ii = 1:size(ff.h_axes,2)
        tR = Rs{an,ii};
        NCells = size(tR.sp_rasters,3);
        tmRR = mRR{an,ii};
        axes(ff.h_axes(1,ii));
        binwidth = tR.bin_width;
        totx = mxsz(ii)*binwidth;
        xs2 = round(totx,0); xs1 = round(max(totx)/2,1);
        xdata = [0 xs1 xs2];
        ydata = [1 size(tmRR,1)];
        imagesc(xdata,ydata,tmRR,[m_mRR M_mRR]); set(gca,'Ydir','normal');
        if ii == 1
            ylabel('Cell #');
        end
        set(gca,'YTick',[],'XTick',[]);

        ylims = ylim;
        if ii == 1
            textstr = sprintf('C%d  %d/%d',conf(ii),floor(ylims(2)),NCells); set_axes_top_text_no_line(ff.hf,gca,textstr,[0.0 -0.1 0.05 0]);
        else
            textstr = sprintf('C%d    %d',conf(ii),floor(ylims(2))); set_axes_top_text_no_line(ff.hf,gca,textstr,[0.0 -0.1 0 0]);
        end
%         textstr = sprintf('C%d',conf(ii)); set_axes_top_text_no_line(ff.hf,gca,textstr,[0 -0.085 0 0]);
        format_axes(gca);
        if ii == 6
          mM = min(tmRR(:)); MM = max(tmRR(:)); 
          [hc,hca] = putColorBar(gca,cbar_p_shift,[m_mRR M_mRR],5,'eastoutside',[0.15 0.18 0.12 0.3]);
        end
        
        axes(ff.h_axes(2,ii));
        tmRR = CRc{an,ii};
        if sum(isnan(tmRR(:))) > 0
            tmRR = fillmissing(tmRR,'linear',2,'EndValues','nearest');
            if sum(isnan(tmRR(:))) > 0
                tmRR = tmRR'; tmRR = fillmissing(tmRR,'linear',2,'EndValues','nearest'); tmRR = tmRR';
            end
        end
        if tim == 0
          if ii <= 3
              totx = 150;
          else
          totx = max(tR.xs);%mxsz(ii)*binwidth;
          end
        else
          if ii > 3
              totx = 15;
          else
          totx = max(tR.xs);%mxsz(ii)*binwidth;
          end
        end
        xs2 = round(totx,0); xs1 = round(max(totx)/2,1);
        xdata = [0 xs2];
        ydata = [1 size(tmRR,1)];
        imagesc(xdata,ydata,tmRR,[m_CRc M_CRc]); set(gca,'Ydir','normal');
        if ii == 1
          if tim
            ylabel('Time (s)');
          else
            ylabel('Dist (cm)');
          end
        end
        set(gca,'YTick',[],'XTick',xdata);

        ylims = ylim;
        textstr = sprintf('%d',ceil(ylims(2)));
%         set_axes_top_text_no_line(ff.hf,gca,textstr,[0 -0.07 0 0]);
        format_axes(gca);
        mM = -0.3;min(tmRR(:)); MM = max(tmRR(:));
        if ii == 6
          [hc,hca] = putColorBar(gca,cbar_p_shift,[m_CRc M_CRc],5,'eastoutside',[0.13 0.18 0.12 0.3]);
        end
        
        axes(ff.h_axes(3,ii));
        tmRR = aCRc{ii};
        if tim
          if ii > 3
              xs1 = 7.5; xs2 = 15;
          else
              totx = mxsz(ii)*binwidth;
              xs2 = round(totx,0); xs1 = round(max(totx)/2,1);
          end
        else
          if ii < 4
              xs1 = 75; xs2 = 150;
          else
              totx = mxsz(ii)*binwidth;
              xs2 = round(totx,0); xs1 = round(max(totx)/2,1);
          end
        end
        xdata = [0 xs1 xs2];
        ydata = [1 size(tmRR,1)];
        imagesc(xdata,ydata,tmRR,[m_aCRc M_aCRc]); set(gca,'Ydir','normal');
        if ii == 1
            if tim
            ylabel('Time (s)');
          else
            ylabel('Dist (cm)');
          end
        end
        set(gca,'YTick',[],'XTick',[0 xs2]);
        if tim
            xlabel('Time (s)');
          else
            xlabel('Dist (cm)');
          end
        
        ylims = ylim;
        textstr = sprintf('%d',ceil(ylims(2)));
%         set_axes_top_text_no_line(ff.hf,gca,textstr,[0 -0.07 0 0]);
        format_axes(gca);
        mM = -0.3;min(tmRR(:)); MM = max(tmRR(:)); 
        if ii == 6
          [hc,hca] = putColorBar(gca,cbar_p_shift,[m_aCRc M_aCRc],5,'eastoutside',[0.13 0.18 0.12 0.3]);
        end
    end
    colormap jet
    ht = set_axes_top_text(gcf,ff.h_axes(1,1),'Air',{0.08,[0.03 -0.05 0 0]});set(ht,'FontWeight','Bold');
    ht = set_axes_top_text(gcf,ff.h_axes(1,2),'No Air',{0.08,[0.03 -0.05 0 0]});set(ht,'FontWeight','Bold');
    ht = set_axes_top_text(gcf,ff.h_axes(1,3),'Air',{0.08,[0.03 -0.05 0 0]});set(ht,'FontWeight','Bold');
    ht = set_axes_top_text(gcf,ff.h_axes(1,4),'No Air',{0.08,[0.03 -0.05 0 0]});set(ht,'FontWeight','Bold');
    save_pdf(ff.hf,mData.pdf_folder,sprintf('pop_vectors_%d_%c.pdf',ntrials,G),600);
