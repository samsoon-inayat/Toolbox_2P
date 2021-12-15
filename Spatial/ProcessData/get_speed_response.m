function ei_C = get_speed_response(ei_C,owr)
for ii = 1:length(ei_C)
    ei_C{ii} = process_one_animal(ei_C{ii},1,owr);
end

function ei = process_one_animal(ei,thr,owr)
% disp(sprintf('%d,%d - %d,%d',length(mOnst),length(mOnse),length(mOffst),length(mOffse)));
binwidths = evalin('base','binwidths');
for pp = 1:length(ei.plane)
    ei.tP = ei.plane{pp}.tP;
    thispFolder = ei.plane{pp}.folder;
    ei.deconv = ei.plane{pp}.tP.deconv;
    
%     rasters = make_rasters(ei,pp,mOnst,mOnse,'time',binwidths);
%     trials = 1:size(rasters.sp_rasters1,1);
%     rasters = findRasterProperties_1(thispFolder,0,'motionOnsets11T',rasters,'time',trials,owr(1:3));
%     ei.plane{pp}.motionOnset_rasters = rasters;
%     
%     rasters = make_rasters(ei,pp,mOffst,mOffse,'time',binwidths);
%     trials = 1:size(rasters.sp_rasters1,1);
%     rasters = findRasterProperties_1(thispFolder,0,'motionOffsets11T',rasters,'time',trials,owr(1:3));
%     ei.plane{pp}.motionOffset_rasters = rasters;
    
%     ei.plane{pp}.acc_response = find_acc_response(ei,pp,owr(5));
    ei.plane{pp}.speed_response = find_speed_response(ei,pp,owr);
end

function out = find_speed_response(ei,pp,owr)
activity = ei.deconv.spSigAll;
speed = ei.b.fSpeed(ei.plane{pp}.b.frames_f);
speed = speed(1:size(activity,2));
out.corr = corr(activity',speed');

min_speed = 1; max_speed = 30;
bin_incr = 1;
bins = min_speed:bin_incr:max_speed;
for ii = 1:(length(bins)-1)
    st = bins(ii);
    se = bins(ii+1);
    bin_centers(ii) = st + bin_incr/2;
    inds = find(speed > st & speed < se);
    cell_act(:,ii) = nanmean(activity(:,inds),2);
end
out.McN = find_cellular_speed_tuning_McN(ei,pp,bins,activity,owr(1));
out.fits = find_cellular_speed_tuning(ei,pp,bin_centers,cell_act,owr(2));
out.FR_vs_speed = cell_act;
out.bin_centers = bin_centers;

function out = find_cellular_speed_tuning_McN(ei,pp,bins,activity,owr)
n = 0;
file_name = fullfile(ei.plane{pp}.folder,sprintf('speed_tuning_McN.mat'));
if exist(file_name,'file') && owr == 0
    out = load(file_name);
    return;
end

speed = ei.b.fSpeed(ei.plane{pp}.b.frames_f);
speed = speed(1:size(activity,2));

hist_speed = hist(speed,bins);
n_hist_s = hist_speed/sum(hist_speed);
% figure(100000);clf;plot(bins,n_hist_s);
r2 = find_cell_speed_tuning(activity,speed,bins,n_hist_s);
nshuffle = 1000;
parfor ii = 1:nshuffle
    shuf_inds = randperm(size(activity,2));
    activity_s = activity(:,shuf_inds);
    r2s(:,ii) = find_cell_speed_tuning(activity_s,speed,bins,n_hist_s);
end
r2_r = repmat(r2,1,nshuffle);
p_vals = sum(r2_r > r2s,2)/nshuffle;
sp_cells = p_vals > 0.99;
inds_nonzero = activity > 0;
speed_tuning_inc = NaN(size(sp_cells,1),length(bins));
speed_tuning_dec = NaN(size(sp_cells,1),length(bins));
speed_cell_resp = double(sp_cells);
for ii = 1:size(sp_cells,1)
    if ~sp_cells(ii)
        continue;
    end
    speeds_when_firing = speed(inds_nonzero(ii,:));
    hist_speed_when_firing = hist(speeds_when_firing,bins);
    n_hist_swf = hist_speed_when_firing/sum(hist_speed_when_firing);
    velocity_tuning = n_hist_swf./n_hist_s;
    [p,S,mu] = polyfit(bins,velocity_tuning,1);
    if p(1) > 0 
        speed_tuning_inc(ii,:) = velocity_tuning;
    end
    if p(1) < 0 
        speed_cell_resp(ii) = -1;
        speed_tuning_dec(ii,:) = velocity_tuning;
    end
%     [f_vt,delta] = polyval(p,bins,S,mu);
%     figure(100000);clf;plot(bins,velocity_tuning,'.');hold on;
%     plot(bins,f_vt);
%     pause;
end
out.speed_tuning_inc = speed_tuning_inc;
out.speed_tuning_dec = speed_tuning_dec;
out.speed_resp = speed_cell_resp;
save(file_name,'-struct','out','-v7.3');
    




function r2 = find_cell_speed_tuning(activity,speed,bins,n_hist_s)
inds_nonzero = activity > 0;
r2 = NaN(size(activity,1),1);
for ii = 1:size(activity,1)
    speeds_when_firing = speed(inds_nonzero(ii,:));
    hist_speed_when_firing = hist(speeds_when_firing,bins);
    n_hist_swf = hist_speed_when_firing/sum(hist_speed_when_firing);
    velocity_tuning = n_hist_swf./n_hist_s;
    [p,S,mu] = polyfit(bins,velocity_tuning,1);
    r2(ii) = 1 - (S.normr/norm(velocity_tuning - mean(velocity_tuning)))^2;
end

function fits = find_cellular_speed_tuning(ei,pp,bcs,cell_act,owr)
n = 0;
max_fr = max(cell_act,[],2);
mean_fr = mean(cell_act,2);

statsetfitnlm = statset('fitnlm');
statsetfitnlm.MaxIter = 1000;
statsetfitnlm.TolFun = 1e-10;
% statsetfitnlm.Display = 'iter';
statsetfitnlm.TolX = statsetfitnlm.TolFun;
statsetfitnlm.UseParallel = 1;
% statsetfitnlm.RobustWgtFun = 'welsch';
% fr = cell_act(mean_fr > 0.1,:);
fr = cell_act;
func_names = {'do_gauss_fit','do_sigmoid_fit','do_linear_fit'};
var_names = {'gauss','sigmoid','linear'};
for ii = 1:length(var_names)
    file_name = fullfile(ei.plane{pp}.folder,sprintf('speed_tuning_%s.mat',var_names{ii}));
    if exist(file_name,'file') && owr == 0
        cmdTxt = sprintf('fits.%s = load(file_name);',var_names{ii});
        eval(cmdTxt);
        continue;
    else
        out = [];
        cmdTxt = sprintf('[out.fitted,~,out.coeffsrs] = %s(bcs,fr,statsetfitnlm,[1 0]);',func_names{ii});
        eval(cmdTxt);
        save(file_name,'-struct','out','-v7.3');
        cmdTxt = sprintf('fits.%s = out;',var_names{ii});
        eval(cmdTxt);
    end
end

% function out = find_acc_response(ei,pp,owr)
% activity = ei.deconv.spSigAll;
% acc = diff(ei.b.speed);
% 
% samplingRate = floor(1/(ei.b.si*1e-6));
% coeffs = ones(1, samplingRate)/samplingRate;
% fAcc = filter(coeffs, 1, acc);
% 
% speed = ei.b.fSpeed(ei.plane{pp}.b.frames_f);
% speed = speed(1:size(activity,2));
% out.corr = corr(activity',speed');
% 
% min_speed = 1; max_speed = 40;
% bin_incr = 1;
% bins = min_speed:bin_incr:max_speed;
% for ii = 1:(length(bins)-1)
%     st = bins(ii);
%     se = bins(ii+1);
%     bin_centers(ii) = st + bin_incr/2;
%     inds = find(speed > st & speed < se);
%     cell_act(:,ii) = nanmean(activity(:,inds),2);
% end
% out.fits = find_cellular_acc_tuning(ei,pp,bin_centers,cell_act,owr);
% out.FR_vs_speed = cell_act;
% out.bin_centers = bin_centers;
% 
% function fits = find_cellular_acc_tuning(ei,pp,bcs,cell_act,owr)
% n = 0;
% max_fr = max(cell_act,[],2);
% mean_fr = mean(cell_act,2);
% 
% statsetfitnlm = statset('fitnlm');
% statsetfitnlm.MaxIter = 1000;
% statsetfitnlm.TolFun = 1e-10;
% % statsetfitnlm.Display = 'iter';
% statsetfitnlm.TolX = statsetfitnlm.TolFun;
% statsetfitnlm.UseParallel = 1;
% % statsetfitnlm.RobustWgtFun = 'welsch';
% % fr = cell_act(mean_fr > 0.1,:);
% fr = cell_act;
% func_names = {'do_gauss_fit','do_sigmoid_fit','do_linear_fit'};
% var_names = {'gauss','sigmoid','linear'};
% for ii = 1:length(var_names)
%     file_name = fullfile(ei.plane{pp}.folder,sprintf('speed_tuning_%s.mat',var_names{ii}));
%     if exist(file_name,'file') && owr == 0
%         cmdTxt = sprintf('fits.%s = load(file_name);',var_names{ii});
%         eval(cmdTxt);
%         continue;
%     else
%         out = [];
%         cmdTxt = sprintf('[out.fitted,~,out.coeffsrs] = %s(bcs,fr,statsetfitnlm,[1 0]);',func_names{ii});
%         eval(cmdTxt);
%         save(file_name,'-struct','out','-v7.3');
%         cmdTxt = sprintf('fits.%s = out;',var_names{ii});
%         eval(cmdTxt);
%     end
% end

