function ei_C = get_speed_response(ei_C,owr)
for ii = 1:length(ei_C)
    ei_C{ii} = process_one_animal(ei_C{ii},1,owr);
end

function ei = process_one_animal(ei,thr,owr)
% disp(sprintf('%d,%d - %d,%d',length(mOnst),length(mOnse),length(mOffst),length(mOffse)));
binwidths = evalin('base','binwidths');
for pp = 1%:length(ei.plane)
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
    ei.plane{pp}.accel_responseZ = find_speed_response(ei,pp,owr);
end

function out = find_speed_response(ei,pp,owr)
raw_activity = ei.deconv.spSigAll;
activity = (raw_activity - mean(raw_activity, 2)) ./ std(raw_activity, [], 2);
speed = ei.b.fSpeed(ei.plane{pp}.b.frames_f);
ts = ei.b.ts(ei.plane{pp}.b.frames_f);
speed = speed(1:size(activity,2));
ts = ts(1:size(activity,2));

accel = diff(speed)./diff(ts);
samplingRate = floor(ei.thorExp.frameRate);
coeffs = ones(1, samplingRate)/samplingRate;
f_accel = filter(coeffs, 1, accel);
activitya = activity(:,1:(end-1));

out.corra = corr(activitya',f_accel');
min_accel = -40; max_accel = 40; bin_incra = 2.5; binsa = min_accel:bin_incra:max_accel;
for ii = 1:(length(binsa)-1)
    st = binsa(ii);
    se = binsa(ii+1);
    bin_centersa(ii) = st + bin_incra/2;
    inds = find(f_accel > st & f_accel < se);
    cell_acta(:,ii) = nanmean(activitya(:,inds),2);
end
out.fits = find_cellular_accel_tuning(ei,pp,bin_centersa,cell_acta,owr(1));
out.FR_vs_accel = cell_acta;
out.bin_centers = bin_centersa;

function fits = find_cellular_accel_tuning(ei,pp,bcs,cell_act,owr)
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
    file_name = fullfile(ei.plane{pp}.folder,sprintf('accel_tuning_%sZ.mat',var_names{ii}));
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

