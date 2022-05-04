function caSig = get_calcium_data(ei,pp)
if ~exist('ei','var')
    ei = evalin('base','ei{3}')
end

n = 0;

fileName = fullfile(ei.plane{pp}.s2p_folder,'ratio_model.mat');
if ~exist(fileName,'file')
    ind = strfind(ei.plane{1}.s2p_folder,'\suite2P');
    fileName = fullfile(ei.plane{pp}.s2p_folder(1:(ind-1)),'ratio_model.mat');
end
matF = matfile(fileName);

caSig = matF.ratio_model';


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