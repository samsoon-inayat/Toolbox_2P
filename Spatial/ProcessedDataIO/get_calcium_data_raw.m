function caSig = get_calcium_data(ei,pp)
if ~exist('ei','var')
    ei = evalin('base','ei{3}')
end

n = 0;

fileName = fullfile(ei.plane{pp}.s2p_folder,'timecourses.mat');
if ~exist(fileName,'file')
    ind = strfind(ei.plane{1}.s2p_folder,'\suite2P');
    fileName = fullfile(ei.plane{pp}.s2p_folder(1:(ind-1)),'timecourses.mat');
end
matF = matfile(fileName);

tcs = matF.tcs;
baselines = repmat(tcs.baseline,size(tcs.raw,1),1);
% caSig = (((tcs.raw-(0.7*tcs.neuropil))-baselines)./baselines)';
caSig = (((tcs.raw-(0*tcs.neuropil))-baselines)./baselines)';


% fileName = fullfile(ei.plane{pp}.s2p_folder,'timecourses.mat');
% if ~exist(fileName,'file')
%     ind = strfind(ei.plane{1}.s2p_folder,'\suite2P');
%     fileName = fullfile(ei.plane{pp}.s2p_folder(1:(ind-1)),'timecourses.mat');
% end
% matF = matfile(fileName);
% 
% tcs = matF.tcs;
% caSig = tcs.raw;


% %%
% fileName = fullfile(ei.plane{pp}.s2p_folder,'Fall.mat');
% if ~exist(fileName,'file')
%     ind = strfind(ei.plane{1}.s2p_folder,'\suite2P');
%     fileName = fullfile(ei.plane{pp}.s2p_folder(1:(ind-1)),'ratio_model.mat');
% end
% matF = matfile(fileName);
% ops = matF.ops;
% %%
% xoff = ops.xoff;
% yoff = ops.yoff;
% figure(100);clf;
% plot(xoff)