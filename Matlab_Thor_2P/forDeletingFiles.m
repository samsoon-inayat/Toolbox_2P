function forDeletingFiles

% There are three or four contexts
% In each context, I want to define air onset and offset, belt markers,
% motion onset and offsets in the intertrial intervals
% so the type of data set I will have will be context definition and then
% in that context I will have different types of markers i.e. air onset and
% offset, belt markers and motion onsets and offsets

ei = evalin('base','ei([1:6])');
startText = 'mutual';
for ii = 1:length(ei)
    tei = ei{ii};
    pFolder = tei.folders.thispFolder;
    dirTxt = sprintf('%s\\%s*.mat',pFolder,startText);
    files = dir(dirTxt);
    for jj = 1:length(files)
        thisFile = files(jj).name;
        fileName = makeName(thisFile,pFolder)
        delete(fileName);
    end
end