function out = myMetrics(x,y,params)

% params = {'no_of_bins_for_MI',10,'no_of_shuffles_for_norm',500,'animal_info',ei,'overwrite_processing',[0 0]};

param = 'no_of_bins_for_MI'; idx = find(strcmp(params,param)); cmdTxt = sprintf('%s = params{idx+1};',param);eval(cmdTxt);
param = 'no_of_shuffles_for_norm'; idx = find(strcmp(params,param)); cmdTxt = sprintf('%s = params{idx+1};',param);eval(cmdTxt);
param = 'animal_info'; idx = find(strcmp(params,param)); cmdTxt = sprintf('%s = params{idx+1};',param);eval(cmdTxt);
param = 'overwrite_processing'; idx = find(strcmp(params,param)); cmdTxt = sprintf('%s = params{idx+1};',param);eval(cmdTxt);
param = 'air_phase'; idx = find(strcmp(params,param)); cmdTxt = sprintf('%s = params{idx+1};',param);eval(cmdTxt);
param = 'configuration'; idx = find(strcmp(params,param)); cmdTxt = sprintf('%s = params{idx+1};',param);eval(cmdTxt);
param = 'variables'; idx = find(strcmp(params,param)); cmdTxt = sprintf('%s = params{idx+1};',param);eval(cmdTxt);
param = 'bin_type'; idx = find(strcmp(params,param)); cmdTxt = sprintf('%s = params{idx+1};',param);eval(cmdTxt);
param = 'trial_type'; idx = find(strcmp(params,param)); cmdTxt = sprintf('%s = params{idx+1};',param);eval(cmdTxt);


file_name = sprintf('metrics_%s_%s_%s_%s_%s.mat',variables,bin_type,trial_type,configuration,air_phase);
folder = fullfile(animal_info.matlab_folder,"unlv"); 
if ~exist(folder,'dir')
    mkdir(folder)
end
file_name = fullfile(folder,file_name);
if overwrite_processing(1) == 1 || ~exist(file_name)
    disp(file_name);
    MI = NaN(size(x,2),3); PC = MI; nobmi = no_of_bins_for_MI; noshffl = no_of_shuffles_for_norm;
    parfor nn = 1:size(x,2)
        [MI(nn,:),PC(nn,:)] = calc_metrics(x(:,nn),y,nobmi,noshffl);
    end
    out.MI = MI; out.PC = PC;
    save(file_name,'-struct','out','-v7.3');
else
    out = load(file_name);
end


% rng(3,'twister');
% [output ~] = info_metrics_S(x, y, no_of_bins_for_MI, [], no_of_shuffles_for_norm);
n = 0;


function [oMI,oPC] = calc_metrics(x,y,nobMI,noshuffles)

xtemp=x(:); xtemp=xtemp(~isnan(xtemp));
y=y(:); y=y(~isnan(x(:)));


%if nshuffles~=0, disp('SHANNON''s Mutual Information...'); end
xbin=NaN(size(xtemp));
%tempedges = quantile(xtemp,linspace(0,1,nobMI+1)); tempedges(end)=tempedges(end)+0.01;
tempedges=[min(xtemp) quantile(xtemp(xtemp>min(xtemp)),linspace(0,1,nobMI))]; 
tempedges(end)=tempedges(end)+0.01;

for b=1:length(tempedges)-1
    xbin(xtemp>=tempedges(b) & xtemp<tempedges(b+1))=b;
end
MI = MutualInformation(y,xbin);
PC = corr(y, xtemp);
if noshuffles == 0
    oMI = MI; oPC = PC;
    return;
else
    aMI = NaN(1,noshuffles);
    aPC = aMI;
    rng(3,'twister');
    for ii = 1:noshuffles
        xtemps = shuffle(xtemp);
        [aMI(ii),aPC(ii)] = calc_metrics(xtemps,y,nobMI,0);
    end
    mu=mean(aMI); si=std(aMI);
    pMI = sum(aMI>MI)/noshuffles;
    zMI =(MI-mu)/si;
    
    mu=mean(aPC); si=std(aPC);
    pPC = sum(aPC>PC)/noshuffles;
    zPC =(PC-mu)/si;

    oMI = [MI zMI pMI]; oPC = [PC zPC pPC];
end