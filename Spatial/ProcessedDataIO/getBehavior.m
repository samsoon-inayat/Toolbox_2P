function tei = getData_py(f,T)
tei = cell(1,1);
for ii = 1:size(T,1)
    try
        f.recordingFolder = cell2mat(T{ii,6});
    catch
        f.recordingFolder = (T{ii,6});
    end
    disp(f.recordingFolder);
    if ~isempty(strfind(f.recordingFolder,'Missing'))
        continue;
    end
    tei{ii}.recordingFolder = f.recordingFolder;
    pp = 1;
    try
        fileName = makeName(sprintf('behavior%d.mat',pp),cell2mat(T{ii,7}));
    catch
        fileName = makeName(sprintf('behavior%d.mat',pp),(T{ii,7}));
    end
    disp(sprintf('Loading 2P plane %d behavior',pp));
    b = load(fileName);
    b = calcBehav(b);
    tei{ii}.b = b;
end
