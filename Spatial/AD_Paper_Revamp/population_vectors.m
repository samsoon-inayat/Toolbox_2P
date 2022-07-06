%% population vector and correlation sensory
while 1
    ff = makeFigureRowsCols(107,[10 3 3.5 2.7],'RowsCols',[3 4],'spaceRowsCols',[0.03 -0.01],'rightUpShifts',[0.1 0.11],'widthHeightAdjustment',[10 -80]);
%     set(gcf,'color','w');    set(gcf,'Position',[10 3 3.5 3.75]);
    MY = 8; ysp = 1; mY = 0; % responsive cells
    stp = 0.15; widths = [0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5]+0.18; gap = 0.16;
    adjust_axes(ff,[mY MY],stp,widths,gap,{'Cell #'});
    
    an = 4;
    o = oC;  
    si = [C1_t_D C2_t_D C3_t_D C4_t_D];
    Rs = o.Rs(:,si);mR = o.mR(:,si);
    ntrials = 70;
    props1 = get_props_Rs(Rs,ntrials);
    good_FR = props1.good_zMI;
    good_FR = good_zMI_MFR_Gauss_C;
    [CRc,aCRc,mRR] = find_population_vector_corr(Rs,mR,good_FR,0);
    for ii = 1:size(ff.h_axes,2)
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
        set(gca,'YTick',[],'XTick',[]);

        ylims = ylim;
        textstr = sprintf('%d',ceil(ylims(2)));
        set_axes_top_text_no_line(ff.hf,gca,textstr,[0 -0.13 0 0]);
        format_axes(gca);
        
        axes(ff.h_axes(2,ii));
        tmRR = CRc{an,ii};
        if ii < 5
            xdata = [0 75 150];
        else
            xdata = [0 7.5 15];
        end
        ydata = [1 size(tmRR,1)];
        imagesc(xdata,ydata,tmRR,[0 1]); set(gca,'Ydir','normal');
        if ii == 1
            ylabel('Distance (cm)');
        end
        set(gca,'YTick',[],'XTick',[]);

        ylims = ylim;
        textstr = sprintf('%d',ceil(ylims(2)));
%         set_axes_top_text_no_line(ff.hf,gca,textstr,[0 -0.07 0 0]);
        format_axes(gca);
        
        axes(ff.h_axes(3,ii));
        tmRR = aCRc{ii};
        if ii < 5
            xdata = [0 75 150];
        else
            xdata = [0 7.5 15];
        end
        ydata = [1 size(tmRR,1)];
        imagesc(xdata,ydata,tmRR,[0 1]); set(gca,'Ydir','normal');
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
    end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('pop_vectors_%d.pdf',ntrials),600);
    
    
    %%
    break;
end


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

