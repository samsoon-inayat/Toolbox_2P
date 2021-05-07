function make_db(T)
for ii = 1:size(T,1)
    rFolder = getT(T{ii,6})
    eip = thorGetExperimentInfo(rFolder);
    nP = getNumberOfPlanes(eip);
    dbFile = fullfile(rFolder,'make_db.m');
    if exist(dbFile,'file')
        continue;
    end
    itf = 1;
    dbFileStr{itf,1} = sprintf('i = 0;');itf = itf+1;
    for kk = 1:nP
        dbFileStr{itf,1} = sprintf('i = i+1;');itf = itf+1;
        dbFileStr{itf,1} = sprintf('db(i).mouse_name    	= ''%d'';',getT(T{ii,1}));itf = itf+1;
        dbFileStr{itf,1} = sprintf('db(i).date                = ''%s'';',datestr(cell2mat(T{ii,2}),'yyyy-mm-dd'));itf = itf+1;
%         dbFileStr{itf,1} = sprintf('db(i).expText             = ''%s'';',thisName);itf = itf+1;
        dbFileStr{itf,1} = sprintf('db(i).expts               = [%d];',kk);itf = itf+1;
        dbFileStr{itf,1} = sprintf('db(i).selectedPlane    	= %d;',kk);itf = itf+1;
        dbFileStr{itf,1} = sprintf('db(i).metaD           = ''%s'';','metaD.xml');itf = itf+1;
        dbFileStr{itf,1} = sprintf('db(i).nplanes       	= 1;');itf = itf+1;
        dbFileStr{itf,1} = sprintf('db(i).BiDiPhase     	= 0;');itf = itf+1;
        dbFileStr{itf,1} = sprintf('db(i).diameter      	= 30;');itf = itf+1;
        dbFileStr{itf,1} = sprintf('db(i).downsamplespace 	= 1;');itf = itf+1;
        dbFileStr{itf,1} = sprintf('db(i).downsampletime 	= 1;');itf = itf+1;
        dbFileStr{itf,1} = sprintf('db(i).downsampletime 	= 1;');itf = itf+1;
%         substr = rFolder((pos+6):end);
%         posls = strfind(substr,'\');
%         pdPaths{ii,kk} = fullfile(pdPath,sprintf('%s\\%d',substr(1:(posls(end)-1)),kk));
    end
    fileID = fopen(dbFile,'w');
    for kk = 1:(itf-1)
        fprintf(fileID,sprintf('%s\n',dbFileStr{kk}));
    end
    fclose(fileID);
end

function t = getT(T)
if iscell(T)
    t = cell2mat(T);
else
    t = T;
end