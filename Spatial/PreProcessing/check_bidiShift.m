function T = check_bidiShift(T,owr)

if ~exist('owr','var')
    owr = 0;
end

for ii = 1:size(T,1)
%     ii
    try
        rFolder  = cell2mat(T{ii,6});
    catch
        rFolder  = (T{ii,6});
    end
    if ~isempty(strfind(rFolder,'Missing'))
        continue;
    end
    rFolder
    eip = thorGetExperimentInfo(rFolder);
    nP = getNumberOfPlanes(eip);
    fileName = 'bidishift.mat';
    try
        pd_path = cell2mat(T{ii,7});
    catch
        pd_path = (T{ii,7});
    end
    fn = fullfile(pd_path,fileName);
    if exist(fn,'file') && owr == 0
        py_pd(ii,1) = 1;
    else
        bidishift = get_bidi_shift(eip,1,100);
        save(fn,'bidishift');
        py_pd(ii,1) = 1;
    end
end

TTemp = table(py_pd,'VariableNames',{'Bi_Di_Shift'});
T = [T TTemp];








% for ii = 1:size(T,1)
%     try
%         rFolder  = cell2mat(T{ii,6});
%     catch
%         rFolder  = (T{ii,6});
%     end
%     eip = thorGetExperimentInfo(rFolder);
%     nP = getNumberOfPlanes(eip);
%     fileName = 'bidishift.mat';
%     bs = [];
%     for pp = 1:nP
%         try
%             pd_path = cell2mat(T{ii,7+pp});
%         catch
%             pd_path = (T{ii,7+pp});
%         end
%         fn = fullfile(pd_path,fileName);
%         if exist(fn,'file')
%             pd(ii,pp) = 1;
%             temp = load(fn);
%             bs(pp) = temp.bidishift;
%         else
%             if isempty(pd_path)
%                 pd(ii,pp) = NaN;
%             else
%                 bidishift = get_bidi_shift(eip,pp,500);
%                 if ~exist(pd_path,'dir')
%                     mkdir(pd_path);
%                 end
%                 save(fn,'bidishift');
%                 bs(pp) = bidishift;
%                 pd(ii,pp) = 1;
%             end
%         end
%     end
%     try
%         pd_path = cell2mat(T{ii,7});
%     catch
%         pd_path = (T{ii,7});
%     end
%     fn = fullfile(pd_path,fileName);
% %     if exist(fn,'file')
% %         py_pd(ii,1) = 1;
% %     else
%         bidishift = bs;
%         if ~exist(pd_path,'dir')
%             mkdir(pd_path);
%         end
%         save(fn,'bidishift');
%         py_pd(ii,1) = 1;
% %     end
% end
% TTemp = table(py_pd,'VariableNames',{'Bi_Di_Shift'});
% T = [T TTemp];
% 
% 
% 
% 
