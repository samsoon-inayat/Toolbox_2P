function T = check_s2p_run(T)

for ii = 1:size(T,1)
    try
        rFolder  = cell2mat(T{ii,6});
    catch
        rFolder  = (T{ii,6});
    end
    if ~isempty(strfind(rFolder,'Missing'))
        continue;
    end
    try
        pd_path = cell2mat(T{ii,7});
    catch
        pd_path = (T{ii,7});
    end
    files = dir(sprintf('%s\\**\\Fall.mat',pd_path));
    if isempty(files)
%         run_s2p(rFolder,pd_path);
        py_s2p(ii,1) = 0;
    else
        py_s2p(ii,1) = 1;
    end
end

TTemp = table(py_s2p,'VariableNames',{'Py_Suite2P'});
T = [T TTemp];


% 
% function [s2p,nm_s2p,py_s2p] = check_s2p_run(T,f)
%     for ii = 1:size(T,1)
%         try
%             rFolder  = cell2mat(T{ii,6});
%         catch
%             rFolder  = (T{ii,6});
%         end
%         eip = thorGetExperimentInfo(rFolder);
%         nP = getNumberOfPlanes(eip);
%         bs = [];
%         for pp = 1:nP
%             try
%                 pd_path = cell2mat(T{ii,7+pp});
%             catch
%                 pd_path = (T{ii,7+pp});
%             end
%             files = dir(sprintf('%s\\F_*plane1.mat',pd_path));
%             if isempty(files)
%                 f.pd_path = pd_path;
%                 try
%                     f.recordingFolder = cell2mat(T{ii,6});
%                 catch
%                     f.recordingFolder = (T{ii,6});
%                 end
% %                 runSuite2P(f,num2str(cell2mat(T{ii,1})),datestr(cell2mat(T{ii,2}),'yyyy-mm-dd'),'',pp,0,NaN);
%                 s2p(ii,pp) = 0;
%             else
%                 s2p(ii,pp) = 1;
%             end
%             files = dir(sprintf('%s\\F_*_proc.mat',pd_path));
%             if isempty(files)
%                 nm_s2p(ii,pp) = 0;
%             else
%                 nm_s2p(ii,pp) = 1;
%             end
% 
%         end
%         try
%             pd_path = cell2mat(T{ii,7});
%         catch
%             pd_path = (T{ii,7});
%         end
%         files = dir(sprintf('%s\\**\\Fall.mat',pd_path));
%         if isempty(files)
%             py_s2p(ii,1) = 0;
%         else
%             py_s2p(ii,1) = 1;
%         end
%     end
% end
