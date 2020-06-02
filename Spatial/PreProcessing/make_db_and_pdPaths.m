function T = make_db_and_pdPaths(T,config,pd_folder_ind)
% pdPath = '\\mohajerani-nas.uleth.ca\storage2\homes\samsoon.inayat\S_Drive\Processed_Data';
% py_pdPath = '\\mohajerani-nas.uleth.ca\storage2\homes\samsoon.inayat\S_Drive\pySuite2p_Processed_Data';
py_pdPath = config.pdFolder{2};
mainCodeFolder = config.mainCodeFolder;
for ii = 1:size(T,1)
%     ii
    if ii == 3
        n = 0;
    end
    try
        rFolder  = cell2mat(T{ii,6});
    catch
        rFolder  = (T{ii,6});
    end
    if ~isempty(strfind(rFolder,'Missing'))
        continue;
    end
    if ~isempty(strfind(cell2mat(T{ii,4}),'10'))
        copyfile(fullfile(mainCodeFolder,'metaD15.xml'),rFolder);
    else
        copyfile(fullfile(mainCodeFolder,'metaD16.xml'),rFolder);
    end
    eip = thorGetExperimentInfo(rFolder);
    nP = getNumberOfPlanes(eip);
    dbFile = fullfile(rFolder,'make_db.m');
%     if exist(dbFile,'file')
%         continue;
%     end
    pos = strfind(rFolder,'\Data\');
    py_pdPaths{ii,1} = fullfile(py_pdPath,rFolder(pos+6:end));
    itf = 1;
    dbFileStr{itf,1} = sprintf('i = 0;');itf = itf+1;
    for kk = 1:nP
        if (T{ii,1}) == 173511
            n = 0;
        end
        dbFileStr{itf,1} = sprintf('i = i+1;');itf = itf+1;
        dbFileStr{itf,1} = sprintf('db(i).mouse_name    	= ''%d'';',(T{ii,1}));itf = itf+1;
        dbFileStr{itf,1} = sprintf('db(i).date                = ''%s'';',datestr(cell2mat(T{ii,2}),'yyyy-mm-dd'));itf = itf+1;
%         dbFileStr{itf,1} = sprintf('db(i).expText             = ''%s'';',thisName);itf = itf+1;
        dbFileStr{itf,1} = sprintf('db(i).expts               = [%d];',kk);itf = itf+1;
        dbFileStr{itf,1} = sprintf('db(i).selectedPlane    	= %d;',kk);itf = itf+1;
        if ~isempty(strfind(cell2mat(T{ii,4}),'10'))
            dbFileStr{itf,1} = sprintf('db(i).metaD           = ''%s'';','metaD15.xml');itf = itf+1;
        else
            dbFileStr{itf,1} = sprintf('db(i).metaD           = ''%s'';','metaD16.xml');itf = itf+1;
        end
        dbFileStr{itf,1} = sprintf('db(i).nplanes       	= 1;');itf = itf+1;
        dbFileStr{itf,1} = sprintf('db(i).BiDiPhase     	= 0;');itf = itf+1;
        dbFileStr{itf,1} = sprintf('db(i).diameter      	= 30;');itf = itf+1;
        dbFileStr{itf,1} = sprintf('db(i).downsamplespace 	= 1;');itf = itf+1;
        dbFileStr{itf,1} = sprintf('db(i).downsampletime 	= 1;');itf = itf+1;
        dbFileStr{itf,1} = sprintf('db(i).downsampletime 	= 1;');itf = itf+1;
        substr = rFolder((pos+6):end);
        posls = strfind(substr,'\');
%         pdPaths{ii,kk} = fullfile(pdPath,sprintf('%s\\%d',substr(1:(posls(end)-1)),kk));
    end
    fileID = fopen(dbFile,'w');
    for kk = 1:(itf-1)
        fprintf(fileID,sprintf('%s\n',dbFileStr{kk}));
    end
    fclose(fileID);
    if ~exist(py_pdPaths{ii,1},'dir')
        mkdir(py_pdPaths{ii,1});
    end
end


TTemp = table(py_pdPaths,'VariableNames',{'Py_ProcessedData'});
T = [T TTemp];
