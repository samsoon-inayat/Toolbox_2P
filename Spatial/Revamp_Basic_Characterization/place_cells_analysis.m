function place_cells_analysis
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
ei = evalin('base','d15'); 
sC = [2 3 4 5 7];
rasterNames = {'air55T','airD','airD','airD','air55T'};
Rs = get_rasters_data(ei,sC,rasterNames);
effective_sampling_rate = 1/0.15;
samplingRate = {'Ca','Ef','Ef','Ef','Ca'};
timeBefore = [5 NaN NaN NaN 5];

%%
for an = 1:5
    cai_sampling_rate = ei{an}.thorExp.frameRate;
    anei = ei{an};
    anRs = Rs(an,:);
    for ii = 1:length(sC)
        sci = sC(ii);
        tRs = anRs{ii};
        thisContext = anei.plane{1}.contexts(sci);
        disp(thisContext.name)
        isCell = logical(tRs.iscell);
        if strcmp(samplingRate{ii},'Ca')
            tempR = tRs.fromFrames.sp_rasters;%(:,:,isCell);
            tempDur = tRs.fromFrames.duration;
            SR = cai_sampling_rate;
        else
            tempR = tRs.sp_rasters1;%(:,:,isCell);
            tempDur = tRs.duration1;
            SR = effective_sampling_rate;
        end
        rasters{an,ii} = find_resp_mdata(tempR,tempDur,timeBefore(ii),SR,thisContext.name,anei,tRs,isCell);
    end
end
[all_ccs,fR] = get_responsive_cells(rasters);
n = 0;
%%
an = 1; cn = [1 2 3 4 5];
plotRasters_compare(200,rasters(an,cn),find(all_ccs{an,3}))
% plotRasters_compare(200,rasters(an,cn),find(all_ccs{an,1}|all_ccs{an,2}|all_ccs{an,3}|all_ccs{an,4}|all_ccs{an,5}))
%%
an = 1; cn = 3;
plotRasters_simple(rasters{an,cn},find(all_ccs{an,cn}),[])
%% visualization of signals
if 1
    meanRs = calc_mean_rasters(rasters(:,[2 3 4]),1:10);
    meanRs = correct_for_size(meanRs);
    
    for an = 1:5
        cell_list = all_ccs{an,4};% & all_ccs{an,3};
        [popVs{an},cpc{an},ccc{an},~] = calc_pop_vector_corr(meanRs(an,:),cell_list,[1 3]);
    end
    graphs = [];
    for an = 1:5
        for cc = 1:3
            graphs{an,cc} = popVs{an}{cc}.popV;
        end
    %     graphs{an,3} = popVsD{an}{1}.popV;
    end

    general_imagesc({1,[1 2 15,8]},graphs)
    return;
end
%%
ssC = [1 2 3];
for an = 1:5
    resp = [];
    for ii = 1:length(ssC)
        tR = rasters{an,ssC(ii)};
        resp(:,ii) = tR.resp.p < 0.05;
    end
    all_resp{an} = resp;
end
%%
average_OI = nan(3,3);
all_OI_mat = [];
for an = 1:5
    pop_map = [];
    OI = [];
    for rr = 1:3
        for cc = 1:3
            rC1 = all_resp{an}(:,rr);
            rC2 = all_resp{an}(:,cc);
            pop_map(rr,cc) = sum(rC1&rC2)/length(rC1);
            OI(rr,cc) = sum(rC1&rC2)/(sum(rC1)+sum(rC2)-sum(rC1&rC2))
        end
    end
    all_pop_map{an} = pop_map;
    mask = triu(ones(size(OI)),1);mask(mask==0) = NaN
    all_OI{an} = OI.*mask;
    if an == 1
        average_OI = all_OI{an};
    else
        average_OI = average_OI + all_OI{an};
    end
    all_OI_mat(:,:,an) = all_OI{an};
end
average_OI = average_OI/5;

general_imagesc({1,[1 2 15 4]},all_OI')
figure(100);clf;imagesc(average_OI); colorbar;
xticks([1 2 3]);yticks([1 2 3])

%%
oi_12 = squeeze(all_OI_mat(1,2,:));
oi_13 = squeeze(all_OI_mat(1,3,:));
oi_23 = squeeze(all_OI_mat(2,3,:));
[h,p,ci,stat] = ttest2(oi_12,oi_13)
[h,p,ci,stat] = ttest2(oi_12,oi_23)
[h,p,ci,stat] = ttest2(oi_13,oi_23)

data = [oi_12 oi_13 oi_23];dataT = array2table(data);dataT.Properties.VariableNames = {'OI12','OI13','OI23'};
within = table([1 2 3]');within.Properties.VariableNames= {'oi_c'};within.oi_c = categorical(within.oi_c);
ra = repeatedMeasuresAnova(dataT,within,0.05);

