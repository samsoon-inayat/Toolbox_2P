function xics = getXicsFromResponse (dfbyf,frames_info)


stimulusFrame = find((frames_info(4):frames_info(5)) == frames_info(1));
srF = stimulusFrame;
erF = srF + 14;
ebF = stimulusFrame - 1;
sbF = ebF - 14;


for jj = 1:size(dfbyf,1)
    thisdfbyf = dfbyf(jj,:);
    thisdfbyf = thisdfbyf - min(thisdfbyf);
    xics.amplitudes(jj) = max(thisdfbyf)-mean(thisdfbyf(sbF:ebF));
    xics.riseTime(jj) = find(thisdfbyf == max(thisdfbyf)) - stimulusFrame;
    xics.areaUnderTheCurve(jj) = trapz(thisdfbyf(srF:erF));
%     xics.areaUnderTheCurve(jj) = getAreaUnderTheCurve(thisdfbyf,stimulusFrame)
end

% function auc = getAreaUnderTheCurve(df,sf)
% 
% n = 0;

