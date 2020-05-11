function findBiDiPhaseMain(f,animal_id,exp_date,brain_location,dbii,ob,prebidishift)

mainDataFolder = f.mainDataFolder; 
tifFolder = f.tifFolder;
mainFolder = f.sDrive;
pmainFolder = f.nsDrive;
pFolder = f.pFolder;

%     animal_id = '183628';
%     exp_date = '2019-06-17';
%     brain_location = '';
try
    mainDataFolder = f.mainDataFolder; 
    db = get_db(mainDataFolder,animal_id,exp_date,'brain_location',brain_location);
catch
    mainDataFolder = f.mainDataFolderList{2}; 
    db = get_db(mainDataFolder,animal_id,exp_date,'brain_location',brain_location);
end
if isempty(db)
    error
end
%     dbii = 1;
%     ob = 1;
%     prebidishift = nan;
if dbii > length(db)
    return;
end
for iii = 1:length(dbii)
    ei{iii} = getData(animal_id,exp_date,dbii(iii),'behavior',0,'twoP',0,'brain_location',brain_location,'overwrite_behavior',ob);
    if ~ei{iii}.zFastEnable & dbii > 1
        disp('Can''t have two planes in this data set');
        return;
    end
    bidishift{iii} = findBiDiPhase(ei{iii},500,prebidishift);
    if isnan(bidishift{iii})
        error;
    end
end

% for iii = 1:length(dbii)
%     if ei{iii}.zFastEnable
%         files = dir(sprintf('%s\\*.tif',ei{iii}.folders.thisTifFolder));
%         if ei{iii}.db.downsampletime < 1
%             if length(files) < ei{iii}.t_timePoints
%                 raw2tif(ei{iii}.rawFile,ei{iii}.folders.thisTifFolder,ei{iii},'frameNumbers',[(1+length(files)):ei{iii}.timePoints],...
%                     'downsamplespace',ei{iii}.db.downsamplespace,'downsampletime',ei{iii}.db.downsampletime,'bidishift',bidishift{iii}(1),...
%                     'plane',ei{iii}.db.selectedPlane);
%             end
%         else
%             if length(files) < ei{iii}.timePoints
%                 raw2tif(ei{iii}.rawFile,ei{iii}.folders.thisTifFolder,ei{iii},'frameNumbers',[(1+length(files)):ei{iii}.timePoints],...
%                     'downsamplespace',ei{iii}.db.downsamplespace,'downsampletime',ei{iii}.db.downsampletime,'bidishift',bidishift{iii}(1),...
%                     'plane',ei{iii}.db.selectedPlane);
%             end
%         end
%     else
%         files = dir(sprintf('%s\\*.tif',ei{iii}.folders.thisTifFolder));
%         if ei{iii}.db.downsampletime < 1
%             if length(files) < ei{iii}.t_timePoints
%                 raw2tif(ei{iii}.rawFile,ei{iii}.folders.thisTifFolder,ei{iii},'frameNumbers',[(1+length(files)):ei{iii}.timePoints],...
%                     'downsamplespace',ei{iii}.db.downsamplespace,'downsampletime',ei{iii}.db.downsampletime,'bidishift',bidishift{iii}(1));
%             end
%         else
%             if length(files) < ei{iii}.timePoints
%                 raw2tif(ei{iii}.rawFile,ei{iii}.folders.thisTifFolder,ei{iii},'frameNumbers',[(1+length(files)):ei{iii}.timePoints],...
%                     'downsamplespace',ei{iii}.db.downsamplespace,'downsampletime',ei{iii}.db.downsampletime,'bidishift',bidishift{iii}(1));
%             end
%         end
%     end
% end
% 
% % exploreData(ei,1,4); % use this function to go through the raw file;
% % rt
% master_file_2017;


function bidishift = findBiDiPhase(ei,nFrames,bidishift)
saveDataFolder = ei.folders.thisTifFolder;
pSaveDataFolder = ei.folders.thispFolder;
fileName = makeName('bidishift.mat',pSaveDataFolder);
if isnan(bidishift)
    if exist(fileName,'file')
        load(fileName);
        return;
    end
    if ei.zFastEnable
        ii = ei.db.selectedPlane;
        getFrames = ii:(ei.zSteps+1):(nFrames*(ei.zSteps+1));
        frames{1} = getFramesFromRaw(ei,getFrames);
    else
        frames{1} = getFramesFromRaw(ei,1:nFrames);
    end

    % framesN = 1:6;
    % frames147 = getFramesFromRaw(ei,framesN);
    % avgFrame = zeros(size(frames(:,:,1)));
    % for ii = 2:3:99
    %     thisFrame = frames(:,:,ii);
    %     avgFrame = avgFrame + thisFrame;
    % end
    % avgFrame = avgFrame/33;
    % figure(103);clf
    % imagesc(avgFrame);
    % axis equal

    for ii = 1%:length(frames)
        aFrame = mean(frames{ii},3);
        bidishift = 0;
        lims = [100 500];
        while 1
            hf = figure(1000);clf;
            set(hf,'Position',get(0,'Screensize')); 
            subplot 121
            imagesc(aFrame);colormap gray;axis equal;
            xlim(lims); ylim(lims);
            title(sprintf('BiDiPhase = %d',bidishift));
            saFrame = ShiftBiDi(bidishift,aFrame,size(aFrame,1),size(aFrame,2));
            subplot 122;
            imagesc(saFrame);colormap gray;axis equal;
            xlim(lims); ylim(lims);
            keyVal = getkey;%pause;
            if keyVal == 115
                if ei.zFastEnable
                    allbds(ii) = bidishift;
                end
                break;
            end
            if keyVal == 29 % Up arrow
                bidishift = bidishift + 1;
            end
            if keyVal == 28 % down arrow
                bidishift = bidishift - 1;
            end
            if keyVal == 27 % down arrow
                pos = get(gcf,'Position');
                pos(3) = pos(3) + 10;
                pos(4) = pos(4) + 10;
                set(gcf,'Position',pos);
            end
            if keyVal == 108 % l
                lims = input('Enter limits ');
                xlim(lims); ylim(lims);
            end
            if keyVal == 113
                error
                break;
            end
        end
    end
    if ei.zFastEnable
        bidishift = allbds;
    end
end
fileName = makeName('bidishift.mat',pSaveDataFolder);
save(fileName,'bidishift');
