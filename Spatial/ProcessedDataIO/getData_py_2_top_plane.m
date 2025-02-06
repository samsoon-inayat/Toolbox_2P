function tei = getData_py_2(T,owrdeconv)

for ii = 1:size(T,1)
    recordingFolder = cell2mat(T{ii,6});
    disp(recordingFolder);
    tei{ii}.recordingFolder = recordingFolder;
    tei{ii}.matlab_folder = cell2mat(T{ii,8});
    tei{ii}.thorExp = thorGetExperimentInfo(recordingFolder);
    nP = getNumberOfPlanes(tei{ii}.thorExp);
    plane{1} = fullfile(cell2mat(T{ii,7}),'suite2P\plane0');
    plane{2} = fullfile(cell2mat(T{ii,7}),'suite2P\plane1');
    matlab_plane{1} = fullfile(cell2mat(T{ii,8}),'suite2P\plane0');
    matlab_plane{2} = fullfile(cell2mat(T{ii,8}),'suite2P\plane1');
    if size(T,2) == 9
        matlab_planeD{1} = fullfile(cell2mat(T{ii,9}),'suite2P\plane0');
        matlab_planeD{2} = fullfile(cell2mat(T{ii,9}),'suite2P\plane1');
    else
        matlab_planeD{1} = fullfile(cell2mat(T{ii,8}),'suite2P\plane0');
        matlab_planeD{2} = fullfile(cell2mat(T{ii,8}),'suite2P\plane1');
    end
    for pp = 1%:length(plane)
        if isempty(plane{pp})
            continue;
        end
        if ~exist(plane{pp},'dir')
            continue;
        end
        sel_plane = plane{pp};
        fileName = makeName(sprintf('behavior%d.mat',pp),cell2mat(T{ii,8}));
        disp(sprintf('Loading 2P plane %d behavior',pp));
        b = load(fileName);
        tei{ii}.plane{pp}.b.frames_f = b.frames_f;
        if pp == 1
            if length(b.frames_f) < tei{ii}.thorExp.timePoints
                tei{ii}.thorExp.timePoints = length(b.frames_f)-1;
                tei{ii}.thorExp.totalFrames = length(b.frames_f)-1;
                tei{ii}.thorExp.streaming_frames = (length(b.frames_f)-1) * (tei{ii}.thorExp.zSteps + 1);
            end
        end
        tei{ii}.plane{pp}.s2p_folder = plane{pp};
        tei{ii}.plane{pp}.folder = matlab_plane{pp};
        tei{ii}.plane{pp}.folderD = matlab_planeD{pp};
        if ~exist(tei{ii}.plane{pp}.folder,'dir')
            mkdir(tei{ii}.plane{pp}.folder);
        end
        recording.s2p_processed_data_folder = cell2mat(T{ii,7}); 
        tP.areCells = get_cluster_iscell(recording,sprintf('plane%d',pp-1));
        disp(sprintf('Loading 2P plane %d spikes',pp));
        deconv = get_cluster_deconv(recording,sprintf('plane%d',pp-1));
        if size(deconv,2) < size(tP.areCells,1)
            tP.areCells = tP.areCells(tP.areCells(:,1)==1,:);
        end
        tP.iscell = tP.areCells;
%         dFbF = get_cluster_timecourses(recording);
        tP.deconv.spSigAll = deconv';
        tei{ii}.plane{pp}.tP = tP;
%         tei{ii}.plane{pp}.ops1{1} = tP.ops;
    end
    try
    b = rmfield(b,'frames_f');
    b = rmfield(b,'frames_r');
    catch
    end
    b = calcBehav(b);
    tei{ii}.b = b;
end
