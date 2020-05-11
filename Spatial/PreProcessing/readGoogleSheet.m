function gS = readGoogleSheet
c = pwd;
p = mfilename('fullpath'); ind = findstr(p,'\');
p = p(1:(ind(end)-1));
cd(p);
system('C:\Users\samsoon.inayat\AppData\Local\Continuum\anaconda3\envs\suite2p\python.exe T:\GitHub\Toolbox_2P\Spatial\PreProcessing\readGoogleSheet.py','-echo')
gSO = load('temp.mat');
cd(c);
gS = convertToLikeExcelSheet(gSO);


function gSO = convertToLikeExcelSheet(gS)
values = gS.values;
for ii = 2:length(values)
    theseValues = values{ii};
    for jj = 1:size(theseValues,1)
        gSO{ii-1,jj} = theseValues(jj,:);
    end
end