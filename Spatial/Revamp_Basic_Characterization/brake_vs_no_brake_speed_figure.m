function brake_vs_no_brake

%% Speed Figure
while 1
    si = [Ar_On Ar_Off Ar_Offc ArL_On ArL_L ArL_Off ArL_Offc Ars_On Ars_Off Ars_Offc];
    Rs = o.Rs(:,si); Rs_MC = o.RsMC(:,si);
    
    for ii = 1:length(ei)
        b1 = ei{ii}.b;
        for jj = 1:10
            alds(ii,jj) = b1.dist(b1.stim_r(jj+10)) - b1.dist(b1.air_puff_r(jj+20));
        end
    end
    ald = round(mean(alds(:)));
    cTxt = rasterNamesTxt(si); 
    cTxt([3 7 10]) = {'3-Arb','4-Arb','5-Arb'};
    % speed
     ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[1 10],'spaceRowsCols',[0.1 0.03],'rightUpShifts',[0.04 0.2],'widthHeightAdjustment',[-32 -400]);
    set(gcf,'color','w'); set(gcf,'Position',[5 5 6.9 0.75]);
%     [Y,E] = discretize(1:49,3);
    all_speeds = []; 
    ind = 1;
    for cn = 1:10
        mean_speed_over_trials = [];
        aThisSpeed = [];
        for an = 1:size(Rs,1)
            thisSpeed = nanmean(Rs{an,cn}.speed);
%             thisSpeed = nanmean(Rs_MC{an,cn}.sp_rasters);
            mean_speed_over_trials(an,:) = thisSpeed;
        end
        mean_speed_over_bins(:,cn) = nanmean(mean_speed_over_trials,2);
            axes(ff.h_axes(1,cn));

        hold on;
        cis = Rs{1,cn}.resp.cis;
        xs = Rs{1,cn}.xs-Rs{1,cn}.xs(cis(1,2)); N = length(xs);
        mspeed = mean(mean_speed_over_trials(:,1:N)); semspeed = std(mean_speed_over_trials(:,1:N))/sqrt(5);
        
        plot(xs,mspeed);
        shadedErrorBar(xs,mspeed,semspeed);
%         changePosition(gca,[0.1 0.15 -0.05 -0.15]);
        if cn == 1 
            put_axes_labels(gca,{'',[0 0 0]},{{'Speed (cm/s)'},[0 5 0]});
        end
        xbTxt = -1; ybTxt = 31;
        text(xbTxt(1),ybTxt+0,cTxt{cn},'FontSize',6);
        ylim([0 25]);
        box off;
        format_axes(gca);
        plot([0 0],[0 30],'m','linewidth',0.25);
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
    %%
    break;
end


%% Speed Figure (12 columns including light arbitrary during conditions 3 and 5)
while 1
    si = [Ar_On Ar_L Ar_Off Ar_Offc ArL_On ArL_L ArL_Off ArL_Offc Ars_On Ars_L Ars_Off Ars_Offc];
    Rs = o.Rs(:,si); Rs_MC = o.RsMC(:,si);
    
    for ii = 1:length(ei)
        b1 = ei{ii}.b;
        for jj = 1:10
            alds(ii,jj) = b1.dist(b1.stim_r(jj+10)) - b1.dist(b1.air_puff_r(jj+20));
        end
    end
    ald = round(mean(alds(:)));
    cTxt = rasterNamesTxt(si); 
%     cTxt([3 7 10]) = {'3-Arb','4-Arb','5-Arb'};
    % speed
     ff = makeFigureRowsCols(107,[1 0.5 4 0.5],'RowsCols',[1 12],'spaceRowsCols',[0.1 0.043],'rightUpShifts',[0.04 0.2],'widthHeightAdjustment',[-45 -400]);
    set(gcf,'color','w'); set(gcf,'Position',[5 5 6.9 0.75]);
%     [Y,E] = discretize(1:49,3);
    all_speeds = []; 
    ind = 1;
    for cn = 1:12
        mean_speed_over_trials = [];
        aThisSpeed = [];
        for an = 1:size(Rs,1)
            thisSpeed = nanmean(Rs{an,cn}.speed);
%             thisSpeed = nanmean(Rs_MC{an,cn}.sp_rasters);
            mean_speed_over_trials(an,:) = thisSpeed;
        end
        mean_speed_over_bins(:,cn) = nanmean(mean_speed_over_trials,2);
            axes(ff.h_axes(1,cn));

        hold on;
        cis = Rs{1,cn}.resp.cis;
        xs = Rs{1,cn}.xs-Rs{1,cn}.xs(cis(1,2)); N = length(xs);
        mspeed = mean(mean_speed_over_trials(:,1:N)); semspeed = std(mean_speed_over_trials(:,1:N))/sqrt(5);
        
        plot(xs,mspeed);
        shadedErrorBar(xs,mspeed,semspeed);
%         changePosition(gca,[0.1 0.15 -0.05 -0.15]);
        if cn == 1 
            put_axes_labels(gca,{'',[0 0 0]},{{'Speed (cm/s)'},[0 5 0]});
        end
        xbTxt = -1; ybTxt = 31;
        text(xbTxt(1),ybTxt+0,cTxt{cn},'FontSize',6);
        ylim([0 25]);
        box off;
        format_axes(gca);
        plot([0 0],[0 30],'m','linewidth',0.25);
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
    %%
    break;
end
