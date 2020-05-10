function place_field_analysis_3
% close all
thisCols_all = {[0 0 0],'b','r',[0 0.7 0.3],[0 0.7 0.3],'r',[0 0.7 0.3],[0 0 0],'b','m','r'};

ei = evalin('base','ei');
owsi = 0;
tei = ei([1 2]);

fileName = makeName('trialsRastersMetaData.mat',tei{1}.folders.thispFolder);
trialsRastersMD = load(fileName);
trials = trialsRastersMD.trials;
% trials = {[],[1:15],[16:31],[32:47]};
trl1 = trials{1}; trl2 = trials{2}; trl3 = trials{3}; trl4 = trials{4};
% trl1 = 1:10; trl3 = 11:26; trl4 = 27:42;
if ~isempty(trl1)
    nbb_air = getRasters_air(tei,trl1);
    temp = find_mutual_information(tei,nbb_air,owsi); nbb_air.SI = temp.zMI;
    nbb_belt = getRasters_belt(tei,trl1);
    temp = find_mutual_information(tei,nbb_belt,owsi); nbb_belt.SI = temp.zMI;
end
nbc_air = getRasters_air(tei,trl2);
temp = find_mutual_information(tei,nbc_air,owsi); nbc_air.SI = temp.zMI;
nbc_belt = getRasters_belt(tei,trl2);
temp = find_mutual_information(tei,nbc_belt,owsi); nbc_belt.SI = temp.zMI;
bc_air_belt = getRasters_air_belt(tei,trl3);
temp = find_mutual_information(tei,bc_air_belt,owsi); bc_air_belt.SI = temp.zMI;
bb_air_belt = getRasters_air_belt(tei,trl4);
temp = find_mutual_information(tei,bb_air_belt,owsi); bb_air_belt.SI = temp.zMI;


% selCells{1} = 1:size(nbc_air.rasters,3);
zScoreTh = 20;
% selCells{1} = find(bb_air_belt.SI > zScoreTh);% & nbb_belt.SI < zScoreTh & bc_air_belt.SI > zScoreTh);
selCells{1} = find(nbc_air.SI < 7 & nbc_belt.SI < 7 & bc_air_belt.SI > zScoreTh)% & nbb_belt.SI < zScoreTh & bc_air_belt.SI > zScoreTh);
length(selCells{1})
n = 0;

varNames = {'nbc_air','nbc_belt','bc_air_belt','bb_air_belt'};
for ii = 1:length(varNames)
    cmdTxt = sprintf('[pcw,pcc] = getPlaceCellProps(%s,selCells{1}); %s.pcw = pcw; %s.pcc = pcc;',varNames{ii},varNames{ii},varNames{ii});
    eval(cmdTxt);
end
% 
% varNames = {'bc_air_belt','bb_air_belt'};
% for ii = 1:length(varNames)
%     cmdTxt = sprintf('pfs = getPlaceCellPropsGauss(%s,selCells); %s.pfs = pfs;',varNames{ii},varNames{ii});
%     eval(cmdTxt);
% end

plotRastersMulti({nbc_air,nbc_belt,bc_air_belt,bb_air_belt},selCells{1},0,0);
return;

