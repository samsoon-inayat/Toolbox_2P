function out = find_cell_speed_tuning(activity,speed)
ei = [];
pp = [];
owr = [0 0];

% activity = ei.deconv.spSigAll;
% speed = ei.b.fSpeed(ei.plane{pp}.b.frames_f);
% speed = speed(1:size(activity,2));
out.corr = corr(activity',speed');

min_speed = 1; max_speed = 40;
bin_incr = 1;
bins = min_speed:bin_incr:max_speed;

min_speed = 1; max_speed = 30;
bin_incr = 1;
new_bins = min_speed:bin_incr:max_speed;
for ii = 1:(length(bins)-1)
    st = bins(ii);
    se = bins(ii+1);
    bin_centers(ii) = st + bin_incr/2;
    inds = find(speed > st & speed < se);
    cell_act(:,ii) = nanmean(activity(:,inds),2);
end
out.fits = find_cellular_speed_tuning(ei,pp,bin_centers,cell_act,owr(2));
out.FR_vs_speed = cell_act;
out.bin_centers = bin_centers;

% for ii = 1:size(cell_act,1)
%     [output ~] = info_metrics_S_onlyMI(cell_act(ii,:), [], 4, [], 0);
%     out.bMI(ii) = output.ShannonMI;
%     % out.mMI(ii) = compute_mutual_information(cell_act(ii,:),[],4);
% end

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
% func_names = {'do_gauss_fit'};
% var_names = {'gauss'};
for ii = 1:length(var_names)
    % file_name = fullfile(ei.plane{pp}.folder,sprintf('speed_tuning_%s.mat',var_names{ii}));
    % if exist(file_name,'file') && owr == 0
    %     cmdTxt = sprintf('fits.%s = load(file_name);',var_names{ii});
    %     eval(cmdTxt);
    %     continue;
    % else
        out = [];
        cmdTxt = sprintf('[out.fitted,~,out.coeffsrs] = %s(bcs,fr,statsetfitnlm,[1 0]);',func_names{ii});
        eval(cmdTxt);
        % save(file_name,'-struct','out','-v7.3');
        cmdTxt = sprintf('fits.%s = out;',var_names{ii});
        eval(cmdTxt);
    % end
end