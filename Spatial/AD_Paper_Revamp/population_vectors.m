%% population vector and correlation sensory
    magfac = mData.magfac;
    ff = makeFigureRowsCols(108,[6 3 3.35 2.5],'RowsCols',[3 4],'spaceRowsCols',[0.03 -0.01],'rightUpShifts',[0.1 0.11],'widthHeightAdjustment',[10 -85]);
%     set(gcf,'color','w');    set(gcf,'Position',[10 3 3.5 3.75]);
    MY = 8; ysp = 1; mY = 0; % responsive cells
    stp = 0.15*magfac; widths = ([0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5]+0.11)*magfac; gap = 0.19*magfac;
    adjust_axes(ff,[mY MY],stp,widths,gap,{'Cell #'});
    
    an = 2;o = oA; G = 'A';  
    an = 4;o = oC; G = 'C';
    si = [C1_t_D C2_t_D C3_t_D C4_t_D];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    good_FR = cell_list_op(props1,{'vals','good_zMI'});
%     good_FR = cell_list_op(props1,{'vals'});
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
        textstr = sprintf('Config. %d (C%d)',ii,ii); set_axes_top_text_no_line(ff.hf,gca,textstr,[0 -0.093 0 0]);
        format_axes_b(gca);
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
        format_axes_b(gca);
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
        format_axes_b(gca);
        mM = -0.3;min(tmRR(:)); MM = max(tmRR(:)); [hc,hca] = putColorBar(gca,cbar_p_shift,[mM MM],5,'eastoutside',[0.09 0.11 0.09 0.16]);
    end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('pop_vectors_%d_%c.pdf',ntrials,G),600);


%% population vector and correlation sensory for time rasters

    magfac = mData.magfac;
    ff = makeFigureRowsCols(108,[8 3 3.35 2.5],'RowsCols',[3 4],'spaceRowsCols',[0.03 -0.01],'rightUpShifts',[0.1 0.11],'widthHeightAdjustment',[10 -85]);
%     set(gcf,'color','w');    set(gcf,'Position',[10 3 3.5 3.75]);
    MY = 8; ysp = 1; mY = 0; % responsive cells
    stp = 0.15*magfac; widths = ([0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5]+0.11)*magfac; gap = 0.19*magfac;
    adjust_axes(ff,[mY MY],stp,widths,gap,{'Cell #'});
    
    an = 2;o = oA; G = 'A';  
%     an = 4;o = oC; G = 'C';
%     si = [C1_t_D C2_t_D C3_t_D C4_t_D];
    si = [C1_i_T C2_i_T C3_i_T C4_i_T];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 50;
    props1 = get_props_Rs(Rs,ntrials);
    good_FR = cell_list_op(props1,{'good_zMI','good_Gauss'});
    good_FR = cell_list_op(props1,{'vals'});
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
        textstr = sprintf('Config. %d (C%d)',ii,ii); set_axes_top_text_no_line(ff.hf,gca,textstr,[0 -0.093 0 0]);
        format_axes_b(gca);
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
        format_axes_b(gca);
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
        format_axes_b(gca);
        mM = -0.3;min(tmRR(:)); MM = max(tmRR(:)); [hc,hca] = putColorBar(gca,cbar_p_shift,[mM MM],5,'eastoutside',[0.09 0.11 0.09 0.16]);
    end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('pop_vectors_%d_%c.pdf',ntrials,G),600);


%% population vector and correlation sensory
while 1
    ff = makeFigureRowsCols(107,[1 0.5 4 1],'RowsCols',[1 8],'spaceRowsCols',[0.01 -0.02],'rightUpShifts',[0.07 0.26],'widthHeightAdjustment',[10 -410]);
    set(gcf,'color','w');    set(gcf,'Position',[10 3 6.9 1.25]);
    MY = 8; ysp = 1; mY = 0; % responsive cells
    stp = 0.11; widths = [0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5]+0.2; gap = 0.16;
    adjust_axes(ff,[mY MY],stp,widths,gap,{'Cell #'});
    
    an = 1;
    o = oC;
    
    si = [C1_t_D C2_t_D C3_t_D C4_t_D C1_i_T C2_i_T C3_i_T C4_i_T];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 70;
    props1 = get_props_Rs(Rs,ntrials);
    good_FR = props1.good_FR_and_tuned;
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,good_FR,0);
    for ii = 1:length(ff.h_axes)
        tR = Rs{an,ii};
        tmRR = mRR{an,ii};
        axes(ff.h_axes(1,ii));
        if ii < 5
            xdata = [0 75 150];
        else
            xdata = [0 7.5 15];
        end
        ydata = [1 size(tmRR,1)];
        imagesc(xdata,ydata,tmRR,[0 1]); set(gca,'Ydir','normal');
        if ii == 1
            ylabel('Cell #');
        end
        set(gca,'YTick',[]);
        if ii < 5
            xlabel('Distance (cm)');
        else
            xlabel('Time (sec)');
        end
        ylims = ylim;
        textstr = sprintf('%d',ceil(ylims(2)));
        set_axes_top_text_no_line(ff.hf,gca,textstr,[0 -0.07 0 0]);
        format_axes(gca);
    end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('pop_vectors_%d.pdf',ntrials),600);
    %%
    break;
end

