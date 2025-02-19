function [ops,stat,iscell] = get_ops(ei,pp)
if ~exist('ei','var')
    ei = evalin('base','ei{3}')
end

n = 0;

fileName = fullfile(ei.plane{pp}.s2p_folder,'Fall.mat');
if ~exist(fileName,'file')
    ind = strfind(ei.plane{1}.s2p_folder,'\suite2P');
    fileName = fullfile(ei.plane{pp}.s2p_folder(1:(ind-1)),'Fall.mat');
end
matF = matfile(fileName);

ops = matF.ops;
stat = matF.stat;
iscell = matF.iscell;