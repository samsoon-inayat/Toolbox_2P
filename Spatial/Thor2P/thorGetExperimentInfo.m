function ei = thorGetExperimentInfo (ei)
ei_flag = 0;
if isstruct(ei) 
    dataFolder = ei.folders.rawDataFolder;
    files = dir(sprintf('%s/*.abf',dataFolder));
    ei.abf_file = makeName(files.name,dataFolder);
    xmlFile = makeName('Experiment.xml',dataFolder);
    ei_flag = 1;
end

if ischar(ei)
    dataFolder = ei;
    xmlFile = makeName('Experiment.xml',dataFolder);
    clear ei;
    files = dir(sprintf('%s/*.abf',dataFolder));
    if length(files) > 0
        ei.abf_file = makeName(files.name,dataFolder);
    else
        ei.abf_file = 'No abf file found';
        disp('No abf file found');
    end
end

xDoc = xmlread(xmlFile);
tt = xDoc.getElementsByTagName('LSM');
ttt = tt.item(0);
temp = char(ttt.getAttribute('frameRate'));ei.frameRate = str2num(temp); %pos = strfind(temp,' ');ei.frameRate = str2num(temp(1:(pos-1)));
temp = char(ttt.getAttribute('pixelX')); ei.pixelX = str2num(temp);
temp = char(ttt.getAttribute('pixelY')); ei.pixelY = str2num(temp);
temp = char(ttt.getAttribute('widthUM')); ei.widthUM = str2num(temp);
temp = char(ttt.getAttribute('heightUM')); ei.heightUM = str2num(temp);

% get the number of frames
tt = xDoc.getElementsByTagName('Timelapse'); % first get the section in xml file ... here e.g. Timelapse
ttt = tt.item(0);
temp = char(ttt.getAttribute('timepoints')); ei.timePoints = str2num(temp);ei.totalFrames = str2num(temp);

tt = xDoc.getElementsByTagName('ExperimentNotes');
ttt = tt.item(0);
temp = char(ttt.getAttribute('text')); ei.ExperimentNotes = temp;
ei.rawFile = makeName('Image_0001_0001.raw',dataFolder);
if ~exist(ei.rawFile)
    ei.rawFile = makeName('Image_001_001.raw',dataFolder);
end
if ei_flag
%     if ei.db.downsampletime < 1
%         ei.t_frameRate = ei.frameRate * ei.db.downsampletime;
%         ei.t_timePoints = ei.timePoints * ei.db.downsampletime;
%         ei.t_totalFrames = ei.t_timePoints;
%     end
%     if ei.db.downsamplespace < 1
%         ei.t_pixelX = ei.pixelX * ei.db.downsamplespace;
%         ei.t_pixelY = ei.pixelY * ei.db.downsamplespace;
%     end
end    


tt = xDoc.getElementsByTagName('Streaming');
ttt = tt.item(0);
temp = char(ttt.getAttribute('zFastEnable'));ei.zFastEnable = str2num(temp); % to see if multiplane imaging has been done
temp = char(ttt.getAttribute('frames'));ei.streaming_frames = str2num(temp); % to see if multiplane imaging has been done

tt = xDoc.getElementsByTagName('ZStage');
ttt = tt.item(0);
temp = char(ttt.getAttribute('steps'));ei.zSteps = str2num(temp); % to see if multiplane imaging has been done

if ei.zFastEnable
    ei.frameRate = ei.frameRate/(ei.zSteps+1);
    n = 0;
end
