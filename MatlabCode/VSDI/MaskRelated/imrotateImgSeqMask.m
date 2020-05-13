function imrotateImgSeqMask(animalNumber)

mainDataFolder = getMainFolder_Navvab;

for an = 1:length(animalNumber)
    dataFolder = makeName(animalNumber{an},mainDataFolder);    
    matFiles = dir(sprintf('%s\\*.mat',dataFolder));
    for ii = 1:length(matFiles)
        temp = matFiles(ii).name;
        if isempty(strfind(temp,'spon'))
            continue;
        end
        fileName = makeName(temp,dataFolder);
        if ~exist(fileName,'file')
            continue;
        end
        load(fileName);
        for idx = 1:size(DF_F0,3)
            im = DF_F0(:,:,idx);
            im = rot90(im,2);
            DF_F0(:,:,idx) = im;            
        end
        save(fileName,'DF_F0','-v7.3')
    end
end

for an = 1:length(animalNumber)
    dataFolder = makeName(animalNumber{an},mainDataFolder);
    peDataFolder = makeName('pEvoked',dataFolder);
    mask.bigMask = getMask(animalNumber(an),1.1);
    mask.mask = getMask(animalNumber(an));
    Folders = dir(peDataFolder);
    for ii = 3:length(Folders)
        folderName = Folders(ii).name;
        folderPath = makeName(folderName,peDataFolder);
        if ~exist(folderPath,'dir')
            continue;
        end
        subFolders = dir(folderPath);
        tt = 0;
        ImgSeqAvg = zeros(128,128,100);
        for jj = 3:length(subFolders)
            subFolderName = subFolders(jj).name;
            subFolderPath = makeName(subFolderName,folderPath);
            fileName = makeName('ImgSeq.mat',subFolderPath);
            if ~exist(fileName,'file')
                continue;
            end
            load(fileName)
            tt = tt+1;
            for idx = 1:size(ImgSeq,3)
                im = ImgSeq(:,:,idx);
                im = rot90(im,2);
                ImgSeq(:,:,idx) = im;
            end
            save(fileName,'ImgSeq','-v7.3')
            ImgSeqAvg = ImgSeqAvg + ImgSeq;
        end
        ImgSeq = ImgSeqAvg/tt;
        fileName = makeName('ImgSeq.mat',folderPath);
        save(fileName,'ImgSeq','-v7.3')
    end
end
