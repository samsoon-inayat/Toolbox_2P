ntrials = 50;
si = [C1_t_D C2_t_D C3_t_D C4_t_D C1_i_T C2_i_T C3_i_T C4_i_T];
%     si = [C1_i_T C2_i_T C3_i_T C4_i_T];
Rs_C = oC.Rs(:,si); Rs_A = oA.Rs(:,si); mRs_C = oC.mR(:,si); mRs_A = oA.mR(:,si);
props_C = get_props_Rs(oC.Rs(:,si),ntrials); props_A = get_props_Rs(oA.Rs(:,si),ntrials);
%% Speed Figure
o = oC; ei = ei_C;
props = props_C;
    Rs = o.Rs(:,si); Rs_MC = o.RsMC(:,si);
    
%     for ii = 1:length(ei)
%         b1 = ei{ii}.b;
%         for jj = 1:10
%             alds(ii,jj) = b1.dist(b1.stim_r(jj+10)) - b1.dist(b1.air_puff_r(jj+20));
%         end
%     end
%     ald = round(mean(alds(:)));
    cTxt = rasterNamesTxt(si); 
%     cTxt([3 7 10]) = {'3-Arb','4-Arb','5-Arb'};
    % speed
     ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[1 8],'spaceRowsCols',[0.1 0.03],'rightUpShifts',[0.04 0.2],'widthHeightAdjustment',[-32 -400]);
    set(gcf,'color','w'); set(gcf,'Position',[5 5 6.9 0.75]);
%     [Y,E] = discretize(1:49,3);
    all_speeds = []; 
    ind = 1;
    for cn = 1:8
        mean_speed_over_trials = [];
        aThisSpeed = [];
        for an = 1:size(Rs,1)
            thisSpeed = props.speed{an,cn};
%             thisSpeed = nanmean(Rs_MC{an,cn}.sp_rasters);
            mean_speed_over_trials(an,:) = thisSpeed;
        end
        mean_speed_over_bins(:,cn) = nanmean(mean_speed_over_trials,2);
            axes(ff.h_axes(1,cn));

        hold on;
        N = 50;
        xs = 1:50;
        mspeed = mean(mean_speed_over_trials(:,1:N)); semspeed = std(mean_speed_over_trials(:,1:N))/sqrt(5);
        
        plot(xs,mspeed);
        shadedErrorBar(xs,mspeed,semspeed);
%         changePosition(gca,[0.1 0.15 -0.05 -0.15]);
        if cn == 1 
            put_axes_labels(gca,{'',[0 0 0]},{{'Speed (cm/s)'},[0 0 0]});
        end
        xbTxt = -1; ybTxt = 31;
        text(xbTxt(1),ybTxt+0,cTxt{cn},'FontSize',6);
        ylim([0 25]);
        box off;
        format_axes(gca);
%         plot([0 0],[0 30],'m','linewidth',0.25);
%         if cn == 3
%             hx = xlabel('Time (sec)');
%         end
        if cn > 1
            set(gca,'YTick',[])
        end
        format_axes(gca);
%         if ismember(cn,[1 2 3])
%             set(gca,'XTick',[])
%         end
    end
    save_pdf(ff.hf,mData.pdf_folder,sprintf('speeds_345'),600);
 
