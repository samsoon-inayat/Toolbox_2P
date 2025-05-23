function trial_to_trial_Analysis

%% load data

dur_cells_T_R = [dur_cells_T(:,1) dur_cells_I(:,1) dur_cells_T(:,2) dur_cells_I(:,2) dur_cells_T(:,3) dur_cells_I(:,3)];

dis_cells_T_R = [dis_cells_T(:,1) dis_cells_I(:,1) dis_cells_T(:,2) dis_cells_I(:,2) dis_cells_T(:,3) dis_cells_I(:,3)];

dur_dis_cells = [dur_cells_T_R dis_cells_T_R];


ntrials = 40;
si = [Ar_t_T Ar_i_T ArL_t_T ArL_i_T Ars_t_T Ars_i_T Ar_t_D Ar_i_D ArL_t_D ArL_i_D Ars_t_D Ars_i_D];
si = [Ar_t_T Ar_i_T ArL_t_T ArL_i_T Ars_t_T Ars_i_T];% Ar_t_D Ar_i_D ArL_t_D ArL_i_D Ars_t_D Ars_i_D];
si = [Ar_t_D Ar_i_T ArL_t_D ArL_i_T Ars_t_D Ars_i_T];% Ar_t_D Ar_i_D ArL_t_D ArL_i_D Ars_t_D Ars_i_D];
si_names = rasterNamesTxt(si);
siG = si; RsG = o.Rs(:,si); propsG = get_props_Rs(RsG,ntrials); respG = propsG.all;
mRsG = calc_mean_rasters(RsG,1:10);
trials = mat2cell([1:10]',ones(size([1:10]')));
[allRsC,allmRsT] = get_trial_Rs(o,si,1:10);
lnsi = length(si);
%     respDT = combine_distance_time_rasters(o.Rs(:,si(1:3)),o.Rs(:,si(4:6)),ntrials);
disp('Done');

% Overlap Indices ImageSC all
avgProps = get_props_Rs(RsG,[50,100]); 
respG = avgProps.vals;
an  = 1:5; eic = 1; sp = 0; intersect_with_global = 0; only_global = 0;
allresp = []; ind = 1;
all_peakL = []; allresp_trials = []; allpeakL_trials = [];
for cn = 1:length(si)
    mRsCT = allmRsT{cn};
    resp = []; peak_locations = [];
    for rr = 1:size(mRsCT,1)
        respTrials = []; pLsTrials = [];
        for cc = 1:size(mRsCT,2)
            this_mat = mRsCT{rr,cc};
            [~,peakL] = max(this_mat,[],2);
%                 size_tmat(rr,cc) = size(this_mat,2);
            resp{rr,cc} = sum(this_mat,2) > 0;
            if cc == 1
                respTrials = resp{rr,cc};
            else
                respTrials = [respTrials resp{rr,cc}];
            end
            if intersect_with_global
                resp{rr,cc} = resp{rr,cc} & respG{rr,cn};
            end
            if only_global
                resp{rr,cc} = respG{rr,cn};
            end
            if sp == 1
                if cn == 1
                    respSe = respSeL{eic}; resp{rr,cc} = resp{rr,cc} & respSe{rr,1};
                end
                if cn == 2
                    respSe = respSeL{eic}; resp{rr,cc} = resp{rr,cc} & respSe{rr,2};
                end
                if cn == 3
                    respSe = respSeL{eic}; resp{rr,cc} = resp{rr,cc} & respSe{rr,3};
                end
                if cn == 4
                    respSe = respSeA{eic}; resp{rr,cc} = resp{rr,cc} & respSe{rr,1};
                end
                if cn == 5
                    respSe = respSeA{eic}; resp{rr,cc} = resp{rr,cc} & respSe{rr,2};
                end
            end
            peakL(~resp{rr,cc}) = NaN;
            peak_locations{rr,cc} = peakL;
            if cc == 1
                pLsTrials = peakL;
            else
                pLsTrials = [pLsTrials peakL];
            end
            if rr == 1
                txl{ind} = sprintf('C%dT%d',cn,cc);
                ind = ind + 1;
            end
        end
        allresp_trials{rr,cn} = respTrials;
        allpeakL_trials{rr,cn} = pLsTrials;
%             oc(rr,cn) = find_cells_based_on_cluster(cell2mat(resp(rr,:)));
    end
    allresp = [allresp resp]; all_peakL = [all_peakL peak_locations];

end
i_allresp = cell_list_op(allresp,[],'not');

allrespOR = cell_list_op(allresp,[],'or',1);
allrespAND = cell_list_op(allresp,[],'and',1);

pallrespOR = 100*exec_fun_on_cell_mat(allrespOR,'sum')./exec_fun_on_cell_mat(allrespOR,'length');
pallrespAND = 100*exec_fun_on_cell_mat(allrespAND,'sum')./exec_fun_on_cell_mat(allrespAND,'length');

[mparOR,semparOR] = findMeanAndStandardError(pallrespOR);

disp('Done');


[OIo,mOI,semOI,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(allresp,0.5,0.05);
%     [OIo,mOI,semOI,OI_mato,p_vals,h_vals,all_CI,mCI,semCI,all_CI_mat,uni] = get_overlap_index(all_resp_T,0.5,0.05);
 disp('Done');

 %% plot conj comp1 and comp2
clear respRV conjV comp1V comp2V

xticklabels = rasterNamesTxt(si);
for an = 1:5
    respRV(an,:) = diag(all_CI_mat(:,:,an));
    conjV(an,:) = diag(all_CI_mat(:,:,an),1);
    comp1V(an,:) = diag(uni(:,:,an),1);
    comp2V(an,:) = diag(uni(:,:,an),-1);
end
respV = respRV;
mrespV = mean(respV); semrespV = std(respV)./sqrt(5);
respTW = reshape(mrespV,10,(size(respRV,2)/10)); mrespAct = respTW'; 
respTW = reshape(semrespV,10,(size(respRV,2)/10)); semrespAct = respTW';

respV = conjV;
mrespV = mean(respV); semrespV = std(respV)./sqrt(5);
mrespV(10:10:(size(respRV,2)-1)) = NaN; mrespV(size(respRV,2)) = NaN; mrespV(isnan(mrespV)) = []; 
semrespV(10:10:(size(respRV,2)-1)) = NaN; semrespV(size(respRV,2)) = NaN; semrespV(isnan(semrespV)) = []; 
respTW = reshape(mrespV,9,(size(respRV,2)/10)); mconjAct = respTW'; 
respTW = reshape(semrespV,9,(size(respRV,2)/10)); semconjAct = respTW';

respV = comp1V;
mrespV = mean(respV); semrespV = std(respV)./sqrt(5);
mrespV(10:10:(size(respRV,2)-1)) = NaN; mrespV(size(respRV,2)) = NaN; mrespV(isnan(mrespV)) = []; 
semrespV(10:10:(size(respRV,2)-1)) = NaN; semrespV(size(respRV,2)) = NaN; semrespV(isnan(semrespV)) = []; 
respTW = reshape(mrespV,9,(size(respRV,2)/10)); mcomp1Act = respTW'; 
respTW = reshape(semrespV,9,(size(respRV,2)/10)); semcomp1Act = respTW';

respV = comp2V;
mrespV = mean(respV); semrespV = std(respV)./sqrt(5);
mrespV(10:10:(size(respRV,2)-1)) = NaN; mrespV(size(respRV,2)) = NaN; mrespV(isnan(mrespV)) = []; 
semrespV(10:10:(size(respRV,2)-1)) = NaN; semrespV(size(respRV,2)) = NaN; semrespV(isnan(semrespV)) = []; 
respTW = reshape(mrespV,9,(size(respRV,2)/10)); mcomp2Act = respTW'; 
respTW = reshape(semrespV,9,(size(respRV,2)/10)); semcomp2Act = respTW';

mrespActL = NaN;  semrespActL = NaN; 
xtl = {''}; trialsStr = cellfun(@num2str,trials','un',0);
xaxL = NaN; xticks = []; xtickL =[];
for ii = 1:(size(respRV,2)/10)
    mrespActL = [mrespActL mrespAct(ii,:) NaN]; semrespActL = [semrespActL semrespAct(ii,:) NaN];
    xaxL = [xaxL 1:10 NaN];
    xtl = [xtl trialsStr {''}];
end
xax = 1:length(mrespActL); 
for ii = 1:length(mrespActL)
    if xaxL(ii) == 1 || xaxL(ii) == 9
        xticks = [xticks xax(ii)];
        xtickL = [xtickL xtl(ii)];
    end
end
rlcolor = [0.75 0.75 0.75];
hf = figure(100);clf;
set(hf,'units','inches','position',[5 5 6.9 1.5]);
%     plot(xax,mrespActL,'color',rlcolor);hold on;
%     plot(xlim,[nanmean(mrespActL) nanmean(mrespActL)],'color',rlcolor);hold on;
plot(xlim,[nanmean(mconjAct(:)) nanmean(mconjAct(:))],'color','m');hold on;
iii=1;
theinds = find(isnan(mrespActL));
for ii = find(isnan(mrespActL))
    plot([ii ii],[4 21],'b-');
    if iii <= (size(respRV,2)/10)
        text(ii+2,21,sprintf('%s',xticklabels{iii}),'FontSize',6);
        indsS = (theinds(iii)+1):(theinds(iii+1)-1);
%             shadedErrorBar(indsS,mrespActL(indsS),semrespActL(indsS),{'color',rlcolor},0.5);
        plot(indsS(1:9),mconjAct(iii,:),'m');
        shadedErrorBar(indsS(1:9),mconjAct(iii,:),semconjAct(iii,:),{'color','m'},0.5);
        plot(indsS(1:9),mcomp1Act(iii,:),'m');
        shadedErrorBar(indsS(1:9),mcomp1Act(iii,:),semcomp1Act(iii,:),{'color','c'},0.5);
        plot(indsS(1:9),mcomp2Act(iii,:),'m');
        shadedErrorBar(indsS(1:9),mcomp2Act(iii,:),semcomp2Act(iii,:),{'color','k'},0.5);
        iii=iii+1;
    end
end
xlim([0 length(mrespActL)+1]); ylim([3 27]);
xlabel('Trial-Pairs');ylabel('Cells (%)');box off;
set(gca,'xtick',xticks,'xticklabel',xtickL);
legs = {'Conjunctive Cells      ','Complementary Cells 1','Complementary Cells 2',[9.5 0.1 25 0.2]}; 
putLegendH(gca,legs,{'m','c','k'},'sigR',{[],'anova',[],6});
format_axes(gca);
changePosition(gca,[-0.08 0.1 0.17 -0.1]);
save_pdf(hf,mData.pdf_folder,sprintf('tria_to_trial_unique.pdf'),600);

%% plot average response percentage
mrespV = mean(respRV); semrespV = std(respRV)./sqrt(5);
respTW = reshape(mrespV,10,(size(respRV,2)/10)); mrespAct = respTW';
respTW = reshape(semrespV,10,(size(respRV,2)/10)); semrespAct = respTW';
mrespActL = NaN;  semrespActL = NaN; 
xtl = {''}; trialsStr = cellfun(@num2str,trials','un',0);
xaxL = NaN; xticks = []; xtickL =[];
for ii = 1:(size(respRV,2)/10)
    mrespActL = [mrespActL mrespAct(ii,:) NaN];
    semrespActL = [semrespActL semrespAct(ii,:) NaN];
    xaxL = [xaxL 1:10 NaN];
    xtl = [xtl trialsStr {''}];
end
xax = 1:length(mrespActL); 
for ii = 1:length(mrespActL)
    if xaxL(ii) == 1 || xaxL(ii) == 10
        xticks = [xticks xax(ii)];
        xtickL = [xtickL xtl(ii)];
    end
end

hf = figure(100);clf;
set(hf,'units','inches','position',[5 5 6.9 1]);
plot(xax,mrespActL,'k');hold on;
plot(xlim,[nanmean(mrespActL) nanmean(mrespActL)],'m');
iii=1;
theinds = find(isnan(mrespActL));
for ii = find(isnan(mrespActL))
    plot([ii ii],[11 26],'b-');
    if iii <= (size(respRV,2)/10)
    text(ii+2,55,sprintf('%s',xticklabels{iii}),'FontSize',6);
    indsS = (theinds(iii)+1):(theinds(iii+1)-1);
    shadedErrorBar(indsS,mrespActL(indsS),semrespActL(indsS),{},0.5)
    iii=iii+1;
    end
end
xlim([0 length(mrespActL)+1]); ylim([10 60]);
xlabel('Trials');ylabel('Cells (%)');box off;
set(gca,'xtick',xticks,'xticklabel',xtickL);
%     legs = {'Responsive Cells',[9.5 0.1 34 0.2]}; 
%     putLegendH(gca,legs,{'k'},'sigR',{[],'anova',[],6});
format_axes(gca);
changePosition(gca,[-0.08 0.1 0.17 -0.1]);
save_pdf(hf,mData.pdf_folder,sprintf('tria_to_trial_unique.pdf'),600);

%% heatmap conj
hf = get_figure(6,[8 3 3.25 3.25]);
%     allresp_sp = [allresp resp_speed(:,4)];
allresp_sp = allresp(1:5,:) ;
[mOI,semOI] = heatmap_conj_comp(gca,allresp_sp,2,{si,rasterNamesTxt});
save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_mean_tu_spatial.pdf'),600);

figdim = 1;
hf = get_figure(7,[5 2 figdim+0.1 figdim]);
minI = min(semOI(:)); maxI = max(semOI(:));
sz = size(mOI,1);
oM = ones(size(mOI));mask1 = (triu(oM,0) & tril(oM,0)); mOI(mask1==1) = NaN;
imAlpha=ones(size(mOI));    %imAlpha(isnan(mask))=0.25; 
imAlpha(mask1 == 1) = 0;
im1 = imagesc(semOI,[minI,maxI]);    im1.AlphaData = imAlpha;
set(gca,'color',0.5*[1 1 1]);    colormap parula;    %axis equal
format_axes(gca);
set_axes_limits(gca,[0.5 sz+0.5],[0.5 sz+0.5]);
set(gca,'xtick',1:length(txl),'ytick',1:length(txl),'xticklabels',[],'yticklabels',[],'Ydir','reverse'); xtickangle(75);
changePosition(gca,[0.0 0 -0.05 0]);
set(gca,'Ydir','normal');
box off
colormap jet
save_pdf(hf,mData.pdf_folder,sprintf('OI_Map_sem.pdf'),600);

%% trial by trial recruitment and position in sequence
    an = 1;
    cn = 1;
    resps = [];
    for trN = 1:120
        resps(:,trN) = allresp{an,trN};
        pks(:,trN) = all_peakL{an,trN};
    end
    sresps = sum(resps,2);
    
    psresps = 100*sresps/size(resps,2);
    
%     [clus_inds,clus_centers] = kmeans(psresps,5);
    
    [sr,sri] = sort(sresps);
    sort_resps = resps(sri,:);
%     [uvs,yvs] = trial_collections(resps);
    %% align by peaks then view
    s_cell_list = [];
    for ii = 1:120
    peak1 = pks(:,ii);
    [sr,sri] = sort(peak1);
    s_cell_list(:,ii) = sri;
    sort_resps = resps(sri,:);
    figure(1000);clf;imagesc(sort_resps)
    pause(0.01);
    end
    figure(1000);clf;imagesc(s_cell_list);colorbar
    %%
    sort_resps = resps(clus_inds == 3,:);
    figure(1000);clf;imagesc(sort_resps)
    %% finding all cells that were recruited in different configurations or their phases
    st = 1:10:(length(si)*10);
    et = 10:10:(length(si)*10);
    for ii = 1:length(si)
        t_resp = resps(:,st(ii):et(ii));
        c_ids = [];
        for jj = 1:10
            tt_resp = t_resp(:,jj);
            c_ids = [c_ids;find(tt_resp)];
        end
        ucids(ii) = length(unique(c_ids));
    end
    
    %%
    pLs = all_peakL{an,cn};
    resp_cells_pLs = pLs(resps);
    
    presp = find_percent(allresp);
    varC = presp;
    [within,dvn,xlabels,awithinD] = make_within_table({'Cond','TI','Tr'},[2,2,10]);
    dataT = make_between_table({varC},dvn);
    ra = RMA(dataT,within,{0.05,{'bonferroni','hsd'}});
    ra.ranova
    print_for_manuscript(ra)
    %% Heatmap of active cells
    st = 1:10:(length(si)*10);
    et = 10:10:(length(si)*10);
    for ii = 1:length(si)
        t_resp = allresp(:,st(ii):et(ii));
        allresp_OR(:,ii) = cell_list_op(t_resp,[],'or',1);
    end
%     si = [
    hf = get_figure(6,[8 3 3.25 3.25]);
    [mOI,semOI] = heatmap_conj_comp(gca,allresp_OR,1,{si,rasterNamesTxt,1});
    
     %% for all active cells in a trial, I want to see the peak locations and average shift of peak locations for animals
    all_var = [];
    out = find_shifts(allresp_trials,allpeakL_trials,[],si_names); % last argument empty means all activated cells in a trial
    all_var = [all_var out.fshiftsF out.fshiftsB out.fshiftsZ];
    
    varC = all_var;
    [within,dvn,xlabels,awithinD] = make_within_table({'TD','PH','Cnds','TI'},[2,3,3,2]);
    dataT = make_between_table({varC},dvn);
    ra = RMA(dataT,within,{0.05,{'bonferroni'}});
    ra.ranova
    print_for_manuscript(ra)
    %%
    for sii = 1:length(si)
        distrib = [];
%         si = 1;
        for an = 1:5
            distrib(an,:) = out.dists{an,sii};
        end
        mdist = mean(distrib);
        semdist = std(distrib)/sqrt(5);
        figure(1000);clf;
        plot(out.xs,mdist);
        errorbar(out.xs,mdist,semdist);
        title(sprintf('%s',si_names{sii}));
        pause;
    end
    %% for all active cells in a trial, I want to see the peak locations and average shift of peak locations for animals
    all_var = [];
    out = find_shifts(allresp_trials,allpeakL_trials,dur_dis_cells);
    all_var = [all_var out.fshiftsF out.fshiftsB out.fshiftsZ];
    
    varC = all_var;
    [within,dvn,xlabels,awithinD] = make_within_table({'TD','PH','Cnds','TI'},[2,3,3,2]);
    dataT = make_between_table({varC},dvn);
    ra = RMA(dataT,within,{0.05,{'bonferroni'}});
    ra.ranova
    print_for_manuscript(ra)
    
    %%
    redF = [2]; redV = {2};
%     redF = [1]; redV = {2};
    [dataTR,withinR] = reduce_within_between(dataT,within,redF,redV);
    ra = RMA(dataTR,withinR,{0.05,{'bonferroni'}});
    ra.ranova
    print_for_manuscript(ra)
   
     %% all cells
    all_var = [];
    out = find_shifts(allresp_trials,allpeakL_trials,dur_dis_cells);
%     all_var = [all_var out.ens.shift_trials_cells];
    all_var = [all_var out.ens.f_shifts];
    
    varC = all_var;
    [within,dvn,xlabels,awithinD] = make_within_table({'TD','PH','Cnds','TI'},[2,3,3,2]);
%     [within,dvn,xlabels,awithinD] = make_within_table({'Cnds','TI'},[3,2]);
    dataT = make_between_table({varC},dvn);
    ra = RMA(dataT,within,{0.05,{'bonferroni'}});
    ra.ranova
    print_for_manuscript(ra)
   
    %%
    redF = [2]; redV = {3};
%     redF = [1]; redV = {2};
    [dataTR,withinR] = reduce_within_between(dataT,within,redF,redV);
    ra = RMA(dataTR,withinR,{0.05,{'bonferroni'}});
    ra.ranova
    print_for_manuscript(ra)
    %%
    an = 5;
    big_mR = [];
    for ii = 1:size(mRsG,2)
        big_mR = [big_mR mRsG{an,ii}];
    end
    
    CRc = corr(big_mR);
    figure(1000);clf;
    subplot(2,1,1);imagesc(big_mR);colorbar;set(gca,'Ydir','Normal');
    subplot(2,1,2);imagesc(CRc);colorbar;set(gca,'Ydir','Normal');