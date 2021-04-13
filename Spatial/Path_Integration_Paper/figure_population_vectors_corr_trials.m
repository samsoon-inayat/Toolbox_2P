function figure1_Distributions
mData = evalin('base','mData'); colors = mData.colors; sigColor = mData.sigColor; axes_font_size = mData.axes_font_size;
filename = fullfile(mData.pd_folder,'pop_corr_mat_trials_AD.mat');
pop_corr = load(filename);
pop_corr_C = pop_corr.all_out_C([1:5]);
pop_corr_A = pop_corr.all_out_A([1:5]);
num_animals_C = size(pop_corr_C{1}.avg_C_conds{1},3);
num_animals_A = size(pop_corr_A{1}.avg_C_conds{1},3);

for ii = 1:length(pop_corr_C)
    this_pop = pop_corr_C{ii};
    avgpc = [];
    max_pos = [];
    for cc = 1:size(this_pop.avg_C_conds,2)
        for an = 1:size(this_pop.avg_C_conds{1},3)
            corrmat = this_pop.avg_C_conds{cc}(:,:,an);
            avgpc(an,cc) = find_average_vals(corrmat);
            mR = this_pop.mean_rasters{an,cc};
            [~,max_pos{an,cc}] = max(mR,[],2);
        end
    end
    avgpc_C_trials(:,:,ii) = avgpc;
    all_max_pos_C{ii} = max_pos;
end

for an = 1:size(all_max_pos_C{1},1)
    tmdiff = [];
    for cc = 1:size(all_max_pos_C{1},2)
        tmp = [];
        for ii = 1:length(all_max_pos_C)
            tmp(:,ii) = all_max_pos_C{ii}{an,cc};
        end
        tmdiff = [tmdiff mean(diff(tmp(:,[1 2]),1,2))];
        tmdiff = [tmdiff mean(diff(tmp(:,[2 5]),1,2))];
    end
    var_C(an,:) = tmdiff;
end

for ii = 1:length(pop_corr_A)
    this_pop = pop_corr_A{ii};
    avgpc = [];
    max_pos = [];
    for cc = 1:size(this_pop.avg_C_conds,2)
        for an = 1:size(this_pop.avg_C_conds{1},3)
            corrmat = this_pop.avg_C_conds{cc}(:,:,an);
            avgpc(an,cc) = find_average_vals(corrmat);
            mR = this_pop.mean_rasters{an,cc};
            [~,max_pos{an,cc}] = max(mR,[],2);
        end
    end
    avgpc_A_trials(:,:,ii) = avgpc;
    all_max_pos_A{ii} = max_pos;
end


for an = 1:size(all_max_pos_A{1},1)
    tmdiff = [];
    for cc = 1:size(all_max_pos_A{1},2)
        tmp = [];
        for ii = 1:length(all_max_pos_A)
            tmp(:,ii) = all_max_pos_A{ii}{an,cc};
        end        
        tmdiff = [tmdiff mean(diff(tmp(:,[1 2]),1,2))];
        tmdiff = [tmdiff mean(diff(tmp(:,[2 5]),1,2))];
    end
    var_A(an,:) = tmdiff;
end

% 
% for an = 1:size(avgpc_C_trials,1)
%     col = 1;
%     for cc = 1:4
%         for tt = 1:length(pop_corr_C)
%             pc_C(an,col) = avgpc_C_trials(an,cc,tt);
%             col = col + 1;
%         end
%     end
% end
% 
% for an = 1:size(avgpc_A_trials,1)
%     col = 1;
%     for cc = 1:4
%         for tt = 1:length(pop_corr_C)
%             pc_A(an,col) = avgpc_A_trials(an,cc,tt);
%             if an == 1
%                 varNames{col} = sprintf('C%dT%d',cc,tt);
%                 temp_tcolors{col} = colors{cc};
%                 xticklabels{col} = sprintf('C%d-TG%d',cc,tt);
%             end
%             col = col + 1;
%         end
%     end
% end

for an = 1:size(avgpc_A_trials,1)
    col = 1;
    for cc = 1:4
        for tt = 1:2%length(pop_corr_C)
            if an == 1
                varNames{col} = sprintf('C%dT%d',cc,tt);
                temp_tcolors{col} = colors{cc};
                xticklabels{col} = sprintf('C%d-TG%d',cc,tt);
            end
            col = col + 1;
        end
    end
end
pc_C = var_C;
pc_A = var_A;
n= 0;
%%
dataT = array2table([[ones(size(pc_C,1),1);2*ones(size(pc_A,1),1)] [pc_C;pc_A]]);
dataT.Properties.VariableNames = {'Group',varNames{:}};
numCols = 2;%length(pop_corr_C);
dataT.Group = categorical(int32(dataT.Group))
colVar1 = [ones(1,numCols) 2*ones(1,numCols) 3*ones(1,numCols) 4*ones(1,numCols)];    colVar2 = [1:numCols 1:numCols 1:numCols 1:numCols];

within = table(colVar1',colVar2');
within.Properties.VariableNames = {'Condition','Trials'};
within.Condition = categorical(within.Condition);
within.Trials = categorical(within.Trials);

ra = repeatedMeasuresAnova(dataT,within);

dataT_C = dataT(1:num_animals_C,2:end);
ra_C = repeatedMeasuresAnova(dataT_C,within);

dataT_A = dataT((num_animals_C+1):end,2:end);
ra_A = repeatedMeasuresAnova(dataT_A,within);

tcolors = [temp_tcolors temp_tcolors];
mVar = ra.est_marginal_means.Mean;
semVar = ra.est_marginal_means.Formula_StdErr;
combs = ra.mcs.combs; p = ra.mcs.p; h = ra.mcs.p<0.05;
sum(h)
number_of_bins = 2;%length(pop_corr_C);
xdata = [1:(number_of_bins*4) ((number_of_bins*4)+2):(((number_of_bins*4)+2)+(number_of_bins*4)-1)]; 
colors = mData.colors;
%     hf = figure(5);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 6.9 2],'color','w');
hf = figure(6);clf;set(gcf,'Units','Inches');set(gcf,'Position',[5 7 3.5 1.5],'color','w');
hold on;
ind = 1
maxY = 10;
[hbs,maxY] = plotBarsWithSigLines(mVar,semVar,combs,[h p],'colors',tcolors,'sigColor','k',...
    'maxY',maxY,'ySpacing',0.06,'sigTestName','','sigLineWidth',0.25,'BaseValue',0,...
    'xdata',xdata,'sigFontSize',7,'sigAsteriskFontSize',10,'barWidth',0.7,'sigLinesStartYFactor',0);
for ii = ((number_of_bins*4)+1):length(hbs)
    set(hbs(ii),'facecolor','none','edgecolor',tcolors{ii});
end
% plot([0.5 11],[-0.5 0.5],'linewidth',1.5)
set(gca,'xlim',[0.25 xdata(end)+0.75],'ylim',[-5 maxY+0.1],'FontSize',6,'FontWeight','Bold','TickDir','out');
xticks = xdata; 

set(gca,'xtick',xticks,'xticklabels',xticklabels);
xtickangle(30)
changePosition(gca,[0 0.03 0.08 0]);
put_axes_labels(gca,{[],[0 0 0]},{{'Correlation just','below diagonal'},[0 0 0]});
save_pdf(hf,mData.pdf_folder,sprintf('corr_trials'),600);

function avgC = find_average_vals(corrmat)
maskmat = ones(size(corrmat));
maskmattril = tril(maskmat,-1) & ~tril(maskmat,-2);
avgC = nanmean(corrmat(maskmattril==1));
% [M,c] = contour(imgaussfilt(corrmat,2));
% avgC = max(c.LevelList)-min(c.LevelList);

% len = size(corrmat,1);
% len1 = round(len/5);
% bins = 1:len1:len;
% maskmat = ones(size(corrmat));
% maskmattril = tril(maskmat,-bins(1)) & ~tril(maskmat,-bins(2));
% corr1 = nanmean(corrmat(maskmattril==1));
% maskmattril = tril(maskmat,-(bins(end)-1)) & ~tril(maskmat,-bins(end));
% corr2 = nanmean(corrmat(maskmattril==1));
% avgC = corr1 - corr2;